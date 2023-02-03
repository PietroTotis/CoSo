import portion as Int
from .util import *


class Configuration(object):
    """
    Represents a target configuration

    Attributes
    ----------
    name : str
        user-defined name
    type: str
        can be sequence/subset/partition/composition
    spec: bool
        each type has an alternative, if spec is true we use that:
        sequence -> permutation
        subset -> multisubset
        composition -> multi-composition
        partition -> partition of any size up to n
    domain: SetFormula
        source set for choices
    size: CSize
        length of sequence/size of subset/number of compositions of partitions
    """

    def __init__(self, name, type, domain, size=None):
        self.name = name
        self.domain = domain
        self.type = str(type)
        self.size = size

    def __repr__(self):
        return str(self)

    def __str__(self):
        str = f"{self.type}"
        str += f" ({self.size}) of entity {self.domain} ({self.name})"
        return str

    def labelled(self):
        return self.type in ["sequence", "permutation", "composition"]

    def lvl1(self):
        return self.type in ["sequence", "permutation", "multisubset", "subset"]

    def lvl2(self):
        return not self.lvl1()

    def with_repetition(self):
        return self.type in ["sequence", "multisubset"]


class Variable:
    """
    Generic Variable of a configuration of level i: a set (domain) for level i=1 configuration or
    a set of sets (lifted set) for level i=2 configurations
    """

    def __init__(self, universe, name):
        self._universe = universe
        self.name = name

    @property
    def universe(self):
        return self._universe

    @universe.setter
    def universe(self, universe):
        self._universe = universe

    def __and__(self, rhs):
        """
        Conjoin constraints with a level 1/2 variable

        Args:
            rhs (Variable): variable to conjoin
        """
        raise NotImplementedError()

    def __contains__(self, val):
        """
        Check if an entity (level 1) or a set (level 2) is included in the valid values
        represented by the variable

        Args:
            val (Level i-1): a level i-1 object
        """
        raise NotImplementedError()

    def __hash__(self):
        """
        Required for dicts indexed by variables
        """
        return hash(self.name)

    def __repr__(self):
        """
        Representation of the variable
        """
        return self.__str__()

    def __str__(self):
        """
        String description
        """
        raise NotImplementedError()

    def copy(self):
        """
        Copy to explicitly create a new configuration for different subproblems
        """
        raise NotImplementedError()


class Constraint:
    """
    Generic constraint on a configuration
    """

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """
        String representation
        """
        raise NotImplementedError()

    def set_universe(self, universe):
        """
        Set universe of the constraint

        Args:
            universe (Domain): universe of the problem
        """
        raise NotImplementedError()


class CCounting(Constraint):
    """
    Attributes
    ----------
    formula : SetFormula/CCounting/CSize
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
        """
        Pretty print admissible values in the interval

        Returns:
            str: a nice description with </>/= operators
        """
        if self.values.atomic:
            lb, ub = interval_closed(self.values, 0, P.inf)
            if lb == ub:
                val = f"= {lb}"
            elif lb == 0:
                val = f"<= {ub}"
            elif ub == P.inf:
                val = f">= {lb}"
            else:
                val = f"in [{lb}, {ub}]"
        else:
            val = f"in {self.values}"
        if isinstance(self.formula, CCounting):
            return f"Nr. ({self.formula}) {val}"
        else:
            return f"Nr. {self.formula} {val}"

    def neg(self):
        interval = Int.closedopen(0, Int.inf) - self.values
        return CCounting(self.formula, interval)

    def complement(self, val, n_rest):
        lb = self.values.lower
        ub = self.values.upper
        if lb == ub:
            return Int.singleton(lb - val)
        else:
            comp = self.values
            if ub == Int.inf:
                comp = comp.replace(lower=lb - val, upper=n_rest, right=Int.CLOSED)
            else:
                comp = comp.replace(upper=ub - val)
            return comp

    def copy(self):
        return CCounting(self.formula, self.values)

    def set_universe(self, universe):
        self.formula.universe = universe

    def update_upper_bound(self, ub):
        if self.values.upper == P.inf:
            self.values = self.values.replace(upper=ub, right=P.CLOSED)


class CPosition(Constraint):
    """
    Attributes
    ----------
    struct : str
        name of the target configuration
    pos : int
        position where the property holds
    dformula: SetFormula
        property/allowed set of elements
    """

    def __init__(self, struct, pos, df):
        self.struct = struct
        self.pos = pos
        self.formula = df

    def __str__(self):
        return f"Position {self.pos}: {self.formula}"

    def set_universe(self, universe):
        self.formula.universe = universe


class CSize(Constraint):
    """
    Container for valid cardinalities of lifted sets

    Attributes
    ----------
    name : str
        a string ignored for now
    values : Int
        interval representing set cardinalities allowed by constraints
    universe : Int
        interval representing all possible set cardinalities
    """

    def __init__(self, name, interval):
        self.name = name
        self.values = interval
        self.universe = Int.closedopen(1, Int.inf)

    def __and__(self, rhs):
        inter = self.values & rhs.values
        return CSize(self.name, inter)

    def __eq__(self, rhs):
        return self.values == rhs.values

    def __contains__(self, rhs):
        if isinstance(rhs, CSize):
            return rhs.values in self.values
        else:  # rhs is int
            return rhs in self.values

    def __iter__(self):
        return Int.iterate(self.values, step=1)

    def __str__(self):
        if self.values.lower == self.values.upper:
            return f"size == {self.values.lower}"
        else:
            return f"size in {self.values}"

    def copy(self):
        return CSize(self.name, self.values)

    def neg(self):
        vals = self.universe.difference(self.values)
        neg = CSize(f"not {self.name}", vals)
        return neg

    def set_universe(self, universe):
        self.universe = universe

    def to_list(self):
        return list([i for i in P.iterate(self.values, step=1)])

    def update_upper_bound(self, ub):
        if self.values.upper == P.inf:
            self.values = self.values.replace(upper=ub, right=P.CLOSED)


# class AggFormula(Constraint):
#     """
#     Attributes
#     ----------
#     set : DomainFOrmula
#         set on which the constraint holds
#     op : function {int} -> int
#         aggregate operator
#     values : Interval
#         admissible values
#     """

#     def __init__(self, set, op, interval):
#         self.set = set
#         self.op = op
#         self.values = interval

#     #     self.fix_bound()

#     # def fix_bound(self):
#     #     if self.op == min and self.values.upper>self.set.elems.upper:
#     #         self.values.replace(upper = self.set.elems.upper)

#     def __repr__(self):
#         if self.op == min:
#             s = "min"
#         elif self.op == max:
#             s = "max"
#         else:
#             s = "sum"
#         return f"{s} of {self.set} in {self.values}"
