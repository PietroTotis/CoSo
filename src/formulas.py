import portion
import operator
import copy

from configuration import Domain
from util import *

class AggFormula(object):
    """
    Attributes
    ----------
    set : DomainFOrmula
        set on which the constraint holds
    op : function {int} -> int
        aggregate operator
    values : Interval
        admissible values
    """

    def __init__(self, set, op, interval):
        self.set = set
        self.op = op
        self.values = interval
    #     self.fix_bound()

    # def fix_bound(self):
    #     if self.op == min and self.values.upper>self.set.elems.upper:
    #         self.values.replace(upper = self.set.elems.upper)

    def __repr__(self):
        if self.op == min:
            s = "min"
        elif self.op == max:
            s = "max"
        else:
            s = "sum"
        return f"{s} of {self.set} in {self.values}"

class CountingFormula(object):
    """
    Attributes
    ----------
    formula : DomainFormula/CountingFormula/SizeFormula
        property to count
    values :
        value interval
    """

    def __init__(self, formula, interval):
        self.formula = formula
        self.values = interval

    def __eq__(self, rhs):
        return self.formula == rhs.formula and self.values == rhs.values

    def __str__(self):
        if self.values.lower == self.values.upper:
            val = f"== {self.values.lower}"
        else:
            val = f"in {self.values}"
        if isinstance(self.formula, DomainFormula):
            return f"Nr. {self.formula} {val}"
        else:
            # return str(self.formula)
            return f"Nr. ({self.formula}) {val}"
        
    def __repr__(self):
        return str(self)

    def neg(self):
        interval = portion.closedopen(0,portion.inf) - self.values
        return CountingFormula(self.formula, interval)
    
    def complement(self, val, n_rest):
        lb = self.values.lower
        ub = self.values.upper
        if lb == ub:
            return portion.singleton(lb-val)
        else:
            comp = self.values
            if ub == portion.inf:
                comp = comp.replace(lower= lb - val, upper=n_rest, right=portion.CLOSED)
            else:
                comp = comp.replace(upper = ub - val)
            return comp
    
    def copy(self):
        return CountingFormula(self.formula, self.values)

class DomainFormula(Domain):
    """
    Represents a set by means of intersections/union/complements of the declared sets

    Attributes
    ----------
    name : str/And/Or/Not
        the description of the set operations to generate the set
    elems : DictInterval
        the corresponding elements
    universe : Domain
        the universe of the problem: each complement is computed w.r.t. it
    """

    def __init__(self, formula, elems, universe):
        self.name = formula
        self.elements = elems
        self.universe = universe
        self.n_elements = None

    @staticmethod
    def is_distinguishable(d1, d2):
        # if e is distinguishable truth value in one domain is d1 and the other d2,
        # if any of the two elements is distinguishable then keep distinguishing
        return d1 or d2 

    def __and__(self, rhs):
        dom = self.elements.domain() & rhs.elements.domain()
        comb = self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)
        elems = comb[dom]
        if elems == self.elements:
            f = self.name
        elif elems == rhs.elements:
            f = rhs.name
        else:
            f = And(self, rhs)
        return DomainFormula(f, elems, self.universe)
    
    def __or__(self, rhs):
        comb = self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)
        if comb == self.elements:
            f = self.name
        elif comb == rhs.elements:
            f = rhs.name
        else:
            f = Or(self, rhs)
        return DomainFormula(f, comb, self.universe)

    def __sub__(self, rhs):
        if self.disjoint(rhs):
            f = self.name
        else:
            f = And(self, Not(rhs))
        comb = self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)
        sub_dom = self.elements.domain() - rhs.elements.domain()
        return DomainFormula(f, comb[sub_dom], self.universe)

    def copy(self):
        return DomainFormula(self.name, self.elements, self.universe)

    def neg(self):
        if isinstance(self.name, Not):
            f = self.name.child
        else:
            f = Not(self)
        els = self.universe.elements.domain() - self.elements.domain()
        dom = self.universe.elements[els]
        return DomainFormula(f, dom, self.universe)

    def indistinguishable_subsets(self, dom_formula=None):
        df = dom_formula if dom_formula is not None else self
        f = df.name
        if isinstance(f, And) or isinstance(f, Or): 
            ind_l = self.indistinguishable_subsets(f.left)
            ind_r = self.indistinguishable_subsets(f.right)
            indist = ind_l.union(ind_r)
        elif isinstance(f, Not) :
            indist = self.indistinguishable_subsets(f.child)
        else: # Domain
            if [False] == list(df.elements.values()):
                indist = {df}
            else:
                indist = set()
        return indist


# class InFormula(object):
#     """
#     A choice formula for subsets

#     Attributes
#     ----------
#     struct: str
#         name of the target configuration
#     entity : DomainFormula
#         property that should belong to the set
#     """

#     def __init__(self, struct, dformula):
#         self.struct = struct
#         self.entity = dformula
        
#     def __repr__(self):
#         return str(self)

#     def __str__(self):
#         return f"{self.entity} is in {self.struct}"

class PosFormula(object):
    """
    Attributes
    ----------
    struct : str
        name of the target configuration
    pos : int
        position where the property holds
    dformula: DomainFormula
        property/allowed set of elements
    """
    def __init__(self, struct, pos, df):
        self.struct = struct
        self.pos = pos
        self.dformula = df

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"Position {self.pos}: {self.dformula}"

class SizeFormula(object):
    """
    Container for valid cardinalities of lifted sets

    Attributes
    ----------
    name : str
        a string ignored for now
    values : portion
        interval representing set cardinalities allowed by constraints
    universe : portion
        interval representing all possible set cardinalities
    """

    def __init__(self, name, interval):
        self.name = name
        self.values = interval
        self.universe = portion.closedopen(1, portion.inf)

    def __and__(self, rhs):
        inter = self.values & rhs.values 
        return SizeFormula(self.name, inter)

    def __eq__(self, rhs):
        return self.values == rhs.values

    def __contains__(self, rhs):
        if isinstance(rhs, SizeFormula):
            return rhs.values in self.values
        else: # rhs is int
            return rhs in self.values

    def __iter__(self):
        return portion.iterate(self.values,step=1)

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.values.lower == self.values.upper:
            return f"size == {self.values.lower}"
        else:
            return f"size in {self.values}"
    
    def copy(self):
        return SizeFormula(self.name, self.values)

    def neg(self):
        vals = self.universe.difference(self.values)
        neg = SizeFormula(f"not {self.name}", vals)
        return neg
    