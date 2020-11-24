import portion
import operator

from structure import Domain
from problog.logic import Term


class CountingFormula(object):
    """
    Attributes
    ----------
    formula : DomainFormula/CountingFormula/SizeFormula
        property to count
    op : str
        one in < > >= =< == \=
    val :
        number to count
    """

    def __init__(self, formula, interval):
        self.formula = formula
        self.values = interval
        # self.op = op
        # self._val = val

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

class DomainFormula(object):
    """
    Represents a set by means of intersections/union/complements of the declared sets

    Attributes
    ----------
    universe : str
        the name of the universe of the problem, as each complement is computed w.r.t. that
    formula : ProbLog Term
        the description of the set operations to generate the set
    domain : Domain
        the corresponding elements
    """

    def __init__(self, universe, formula, domain):
        self.universe = universe
        self.formula = formula
        self.domain = domain

    def __and__(self, rhs):
        dom = self.domain & rhs.domain
        if dom == self.domain:
            int_term = self.formula
        elif dom == rhs.domain:
            int_term = rhs.formula
        else:
            int_term = Term("inter", self.formula, rhs.formula)
        return DomainFormula(self.universe, int_term, dom)
    
    def __contains__(self, other):
        return other.domain in self.domain

    def __eq__(self, rhs):
        return self.domain == rhs.domain

    def __hash__(self):
        return hash(self.formula)

    def __or__(self, rhs):
        dom = self.domain | rhs.domain
        if dom in self.domain:
            union_term = self.formula
        elif dom in rhs.domain:
            union_term = rhs.formula
        else:
            union_term = Term("union",  self.formula, rhs.formula)
        return DomainFormula(self.universe, union_term, dom)
        
    def __repr__(self):
        return str(self)

    def __sub__(self, rhs):
        if self.disjoint(rhs):
            sub_formula = self.formula
        else:
            sub_formula = Term("inter", self.formula, Term("not", rhs.formula))
        sub_dom = self.domain - rhs.domain
        return DomainFormula(self.universe, sub_formula, sub_dom)

    def __str__(self):
        if self.domain.size() > 0 :
            str = f"set of {self.domain.name} ({self.domain})"
        else:
            str = f"{self.to_str(self.formula)} (empty)"
        return str

    def copy(self):
        return DomainFormula(self.universe,self.formula, self.domain)

    def disjoint(self, rhs):
        return self.domain.disjoint(rhs.domain)

    def neg(self):
        if self.formula.functor == "not":
            not_term = self.formula.args[0]
        else:
            not_term = Term("not", self.formula)
        dom = self.universe - self.domain
        return DomainFormula(self.universe, not_term, dom)

    def to_str(self, f):
        if f.functor == "inter":
            return " ∧ ".join(map(self.to_str, f.args))
        elif f.functor == "union":
            return " ∨ ".join(map(self.to_str, f.args))
        elif f.functor == "not":
            return f"¬({self.to_str(f.args[0])})"
        else:
            return str(f)       

    def take(self, n):
        c = self.copy()
        c.domain = self.domain.take(n)
        return c

    def size(self):
        return self.domain.size()

class InFormula(object):
    """
    A choice formula for subsets

    Attributes
    ----------
    struct: str
        name of the target structure
    entity : DomainFormula
        property that should belong to the set
    """

    def __init__(self, struct, dformula):
        self.struct = struct
        self.entity = dformula
        
    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.entity} is in {self.struct}"

class PosFormula(object):
    """
    Attributes
    ----------
    struct : str
        name of the target structure
    pos : int
        position where the property holds
    dformula: DomainFormula
        property/allowed set of elements
    """
    def __init__(self, struct, pos, df):
        self.struct = struct
        self.pos = pos.compute_value()
        self.dformula = df

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"Position {self.pos}: {self.dformula}"

class SizeFormula(object):
    """
    Container for allowed cardinalities of lifted sets

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
    
    def neg(self):
        return self.universe.difference(self.values)
    