import portion
import operator

from structure import Domain
from problog.logic import Term


class CountingFormula(object):
    """
    Attributes
    ----------
    dformula : DomainFormula
        property to count
    op : str
        one in < > >= =< == \=
    val :
        number to count
    """

    def __init__(self, formula, op, val):
        self.dformula = formula
        self.op = op
        self._val = val

    def __str__(self):
        return f"Nr. {self.dformula} {self.op} {self._val}"

    def get_operator(self):
        if self.op == ">":
            return operator.gt
        elif self.op == "<":
            return operator.lt
        elif self.op == "=<":
            return operator.le
        elif self.op == ">=":
            return operator.ge
        elif self.op == "==":
            return operator.eq
        elif self.op == "\=":
            return operator.ne

    def invert(self, tot):
        neg_f = self.dformula.neg()
        if self.op == ">":
            return CountingFormula(neg_f, "<", tot - self._val)
        elif self.op == "<":
            return CountingFormula(neg_f, ">", tot - self._val)
        elif self.op == "=<":
            return CountingFormula(neg_f, ">=", tot - self._val)
        elif self.op == ">=":
            return CountingFormula(neg_f, "=<", tot - self._val)
        elif self.op == "==":
            return CountingFormula(neg_f, "\=", tot - self._val)
        else: # self.op == "\=":
            return CountingFormula(neg_f, "==", tot - self._val)
    
    def num(self):
        """
        In the algorithms use <= or >= so adjust the number accordingly
        """
        if self.op == ">":
            return self._val +1
        elif self.op == "<":
            return self._val -1 
        else:
            return self._val

    def update(self, val):
        return CountingFormula(self.dformula, self.op, val)


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

    def __sub__(self, rhs):
        if self.disjoint(rhs):
            sub_formula = self.formula
        else:
            sub_formula = Term("inter", self.formula, Term("not", rhs.formula))
        sub_dom = self.domain - rhs.domain
        return DomainFormula(self.universe, sub_formula, sub_dom)

    def __str__(self):
        str = f"{self.to_str(self.formula)} ({self.domain})"
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

    def __str__(self):
        return f"{self.entity} is in {self.struct}"

class PosFormula(object):
    """
    Attributes:
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

    def __str__(self):
        return f"Position {self.pos}: {self.dformula}"