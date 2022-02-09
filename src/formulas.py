import portion

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

    def __init__(self, formula, elems, universe, name = None):
        if name is None:
            self.name = formula
        else:
            self.name = name
        self.formula = formula
        self.elements = elems
        self.universe = universe
        self.n_elements = None

    def __add__(self, rhs):
        return self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)

    def __and__(self, rhs):
        dom = self.elements.domain() & rhs.elements.domain()
        # comb = self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)
        comb = combine(self, rhs)
        elems = comb[dom]
        if elems == self.elements:
            f = self.formula
        elif elems == rhs.elements:
            f = rhs.formula
        else:
            f = And(self, rhs)
        if self.name == "universe":
            name = rhs.name
        elif rhs.name == "universe":
            name = self.name
        elif self.name == rhs.name:
            name = self.name
        else:
            name = And(str(self), str(rhs))
        return DomainFormula(f, elems, self.universe, name)
    
    def __or__(self, rhs):
        comb = combine(self, rhs)
        if comb.domain() in self.elements.domain():
            f = self.formula
        elif comb.domain() in rhs.elements.domain():
            f = rhs.formula
        else:
            f = Or(self, rhs)
        if self.name == "universe" or rhs.name == "universe":
            name = "universe"
        elif self.name == rhs.name:
            name = self.name
        else:
            name = Or(str(self), str(rhs))
        return DomainFormula(f, comb, self.universe, name)

    def __sub__(self, rhs):
        if self.disjoint(rhs):
            f = self.formula
        else:
            f = And(self, Not(rhs))
        comb = self.elements.combine(rhs.elements, how=DomainFormula.is_distinguishable)
        sub_dom = self.elements.domain() - rhs.elements.domain()
        return DomainFormula(f, comb[sub_dom], self.universe, self.name)

    def copy(self):
        return DomainFormula(self.formula, self.elements, self.universe, self.name)

    def neg(self):
        if isinstance(self.formula, Not):
            f = self.formula.child
        else:
            f = Not(self)
        els = self.universe.elements.domain() - self.elements.domain()
        dom = self.universe.elements[els]
        name = Not(str(self))
        return DomainFormula(f, dom, self.universe, name)

    def indistinguishable_subsets(self, dom_formula=None):
        df = dom_formula if dom_formula is not None else self
        f = df.formula
        if isinstance(f, And) or isinstance(f, Or): 
            ind_l = self.indistinguishable_subsets(f.left)
            ind_r = self.indistinguishable_subsets(f.right)
            indist = ind_l.union(ind_r)
        elif isinstance(f, Not) :
            indist = self.indistinguishable_subsets(f.child)
        else: # Domain
            df = df & self.universe # update indistinguishability from universe
            if [False] == list(df.elements.values()):
                indist = {df}
            else:
                indist = set()
            # if df.name == "universe":

            # else:
            #     df = df & self.universe # update indistinguishability from universe
            #     if [False] == list(df.elements.values()):
            #         indist = {df}
            #     else:
            #         indist = set()
        return indist

    def set_universe(self, universe):
        self.universe = universe
        if isinstance(self.formula, And) or isinstance(self.formula, Or): 
            if isinstance(self.formula.left, DomainFormula):
                self.formula.left.set_universe(universe)
            if isinstance(self.formula.right, DomainFormula):
                self.formula.left.set_universe(universe)
        elif isinstance(self.formula, Not):
            if isinstance(self.formula.child, DomainFormula):
                self.formula.child.set_universe(universe)
        else:
            pass


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
    