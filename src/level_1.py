import portion as P
from configuration import Variable
from util import *


class Multiset(object):
    """
    Represents a domain as a portion (set of intervals) with extra info, i.e. distinguishability, size
    and take operation

    Attributes
    ----------
    name : str
        the label associated to the multiset
    formula: str
        the set-formula describing the multiset (usually same as name)
    elements : IntervalDict
        Intervals corresponding to the entities (mapped to integers) of the multiset.
        Mapped as dict to True if distinguishable False otherwise
    """

    def __init__(self, name, elems):
        self.name = name
        self.formula = name
        self.elements = elems
        self.n_elements = None

    def __contains__(self, val):
        return val.elements.domain() in self.elements.domain()

    def __eq__(self, rhs):
        return self.elements == rhs.elements

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name is not None:
            return str(self.name)
        else:
            s = f"{self.formula}"
            if self.all_indistinguishable() and len(self.elements.domain()) == 1:
                n_elems = sum([1 for _ in P.iterate(self.elements.domain(), step=1)])
                i = self.elements.domain().lower
                e_lab = self.labels.get(i, i)
                s = f"{n_elems}x " + str(e_lab)
        # ok for debug, not for log
        # if self.size() > 0 :
        #     s += f"({list(self.elements.keys())})=|{self.n_elements}|"
        # else:
        #     s += f"(none)"
        return s

    def all_indistinguishable(self):
        return not (True in self.elements.values())

    def max(self):
        return self.elements.domain().upper

    def min(self):
        return self.elements.domain().lower

    def disjoint(self, rhs):
        inter = self & rhs
        return inter.elements.domain().empty

    def size(self):
        if self.n_elements is None:
            s = 0
            for e in self.elements.domain():
                if e.upper == P.inf:
                    return -1
                if not e.empty:
                    if e.left == P.CLOSED and e.right == P.CLOSED:
                        s += e.upper - e.lower + 1
                    elif e.left == P.OPEN and e.right == P.OPEN:
                        s += e.upper - e.lower - 1
                    else:
                        s += e.upper - e.lower
            self.n_elements = s
        return self.n_elements

    def take(self, n):
        """
        returns a subset of itself of size n
        if n > size returns itself
        since we use it for shattering alldiff, assumes that all elements are distinguishable
        """
        iter = P.iterate(self.elements.domain(), step=1)
        subset = P.empty()
        i = 0
        hasNext = True
        while hasNext and i < n:
            try:
                elem = P.singleton(next(iter))
                subset = subset | elem
            except StopIteration:
                hasNext = False
            i += 1
        dist = P.IntervalDict()
        dist[subset] = True
        taken = Multiset(f"{n}x {self.name}", dist)
        return taken


class SetFormula(Multiset, Variable):
    """
    Represents a set by means of intersections/union/complements of the declared sets

    Attributes
    ----------
    name : str/And/Or/Not
        the description of the set operations to generate the set
    elems : DictInterval
        the corresponding elements
    universe : Multiset
        the universe of the problem: each complement is computed w.r.t. it
    """

    def __init__(self, formula, elems, universe, name=None, labels={}):
        Variable.__init__(self, universe, name)
        Multiset.__init__(self, formula, elems)
        if name is None:
            self.name = formula
        else:
            self.name = name
        self.set_labels(labels)

    def __add__(self, rhs):
        return self.elements.combine(rhs.elements, how=is_distinguishable)

    def __and__(self, rhs):
        """
        Intersection of SetFormulas self /\ rhs

        Args:
            rhs (SetFormulas): another SetFormulas

        Returns:
            SetFormulas: SetFormulas representing the intersection
        """
        dom = self.elements.domain() & rhs.elements.domain()
        # comb = self.elements.combine(rhs.elements, how=SetFormula.is_distinguishable)
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
            name = self.simplify_name(elems, rhs, And)
        return SetFormula(f, elems, self.universe, name, self.labels)

    def __or__(self, rhs):
        """
        Disjunction of SetFormulas

        Args:
            rhs (SetFormulas): other SetFormulas

        Returns:
            SetFormulas: union of self with rhs
        """
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
        return SetFormula(f, comb, self.universe, name, self.labels)

    def __sub__(self, rhs):
        """
        Remove objects in rhs from self

        Args:
            rhs (SetFormulas): objects to remove

        Returns:
            SetFormulas: self-rhs
        """
        if self.disjoint(rhs):
            f = self.formula
        else:
            f = And(self, Not(rhs))
        comb = self.elements.combine(rhs.elements, how=is_distinguishable)
        sub_dom = self.elements.domain() - rhs.elements.domain()
        return SetFormula(f, comb[sub_dom], self.universe, self.name, self.labels)

    def copy(self):
        return SetFormula(
            self.formula, self.elements, self.universe, self.name, self.labels
        )

    def neg(self):
        """
        Complement the set w.r.t. universe

        Returns:
            SetFormulas: universe-self
        """
        if isinstance(self.formula, Not):
            f = self.formula.child
            name = str(self)[1:]
        else:
            f = Not(self)
            name = Not(str(self))
        els = self.universe.elements.domain() - self.elements.domain()
        dom = self.universe.elements[els]
        return SetFormula(f, dom, self.universe, name, self.labels)

    def indistinguishable_subsets(self, set_formula=None):
        """
        Collect all partitions of indistinguishable objects contained in the set formula

        Args:
            set_formula (SetFormula, optional): set formula . Defaults to None.

        Returns:
            [portion]: intervals corresponding to the ids of indistinguishable objects
        """
        df = set_formula if set_formula is not None else self
        f = df.formula
        if isinstance(f, And) or isinstance(f, Or):
            ind_l = self.indistinguishable_subsets(f.left)
            ind_r = self.indistinguishable_subsets(f.right)
            indist = ind_l.union(ind_r)
        elif isinstance(f, Not):
            indist = self.indistinguishable_subsets(f.child)
        else:  # Domain
            df = df & self.universe  # update indistinguishability from universe
            vals = list(df.elements.values())
            if [False] == vals:
                indist = {df}
            elif False in vals:
                indist_subsets = set()
                i = 1
                for dom in df.elements.find(False):
                    sub = SetFormula(f"{df.formula}_{i}", df.elements[dom], df.universe)
                    indist_subsets.add(sub)
                    i += 1
                return indist_subsets
            else:
                indist = set()
        return indist

    def set_universe(self, universe):
        """
        Propagate universe info to subsets

        Args:
            universe (Domain): universe
        """
        self.universe = universe
        if isinstance(self.formula, And) or isinstance(self.formula, Or):
            if isinstance(self.formula.left, SetFormula):
                self.formula.left.set_universe(universe)
            if isinstance(self.formula.right, SetFormula):
                self.formula.right.set_universe(universe)
        elif isinstance(self.formula, Not):
            if isinstance(self.formula.child, SetFormula):
                self.formula.child.set_universe(universe)
        else:
            pass

    def set_labels(self, labels):
        """
        Propagate label info to subsets

        Args:
            labels ({int:string}): maps objects ids to label
        """
        self.labels = labels
        if isinstance(self.formula, And) or isinstance(self.formula, Or):
            if isinstance(self.formula.left, SetFormula):
                self.formula.left.set_labels(labels)
            if isinstance(self.formula.right, SetFormula):
                self.formula.right.set_labels(labels)
        elif isinstance(self.formula, Not):
            if isinstance(self.formula.child, SetFormula):
                self.formula.child.set_labels(labels)
        else:
            pass

    def simplify_name(self, elems, rhs, op):
        """
        If a SetFormula resulting from And/Or represents a single element
        then simplify the description using the name of the element

        Args:
            elems (portion): either a single distinguishable element or n indistinguishable copies
            rhs (SetFormula): the other operand
            op (And/Or): the operator to combine the names if simplification is not possible

        Returns:
            str: simplified description
        """
        if len(elems.domain()) == 1:
            if not (True in elems.values()):
                n_elems = sum([1 for _ in P.iterate(elems.domain(), step=1)])
                i = elems.domain().lower
                e_lab = self.labels.get(i, i)
                name = f"{n_elems}x " + str(e_lab)
            else:
                i = elems.domain().lower
                e_lab = self.labels.get(i, i)
                name = f"1x " + str(e_lab)
        else:
            name = op(str(self), str(rhs))
        return name
