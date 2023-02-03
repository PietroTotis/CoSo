import portion as P

from .configuration import Variable
from .util import *


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
        return val.elements.keys() <= self.elements.keys()

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
                # e_lab = self.labels.get(i, i)
                s = f"{n_elems}x " + self.name  # str(e_lab)
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

    # TODO @property
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

    def __init__(self, formula, elems, universe, name=None):
        Variable.__init__(self, universe, name)
        Multiset.__init__(self, formula, elems)
        if name is None:
            self.name = self.universe.get_label(support(elems), str(formula))
        else:
            self.name = name

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
            return self
        elif elems == rhs.elements:
            return rhs
        else:
            f = And(self, rhs)
            default = simplify_name(elems, str(f), self.universe)
            name = self.universe.get_label(dom, default)
            return SetFormula(f, elems, self.universe, name)

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
            return self
        elif comb.domain() in rhs.elements.domain():
            return rhs
        else:
            f = Or(self, rhs)
            default = simplify_name(comb, str(f), self.universe)
            name = self.universe.get_label(support(comb), default)
            return SetFormula(f, comb, self.universe, name)

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
        return SetFormula(f, comb[sub_dom], self.universe, self.name)

    def __iter__(self):
        for e in P.iterate(support(self.elements), step=1):
            yield e

    def copy(self):
        return SetFormula(self.formula, self.elements, self.universe, self.name)

    def get_label(self, key, default=None):
        return self.universe.get_label(key, default)

    def empty(self):
        return support(self.elements).empty

    def neg(self):
        """
        Complement the set w.r.t. universe

        Returns:
            SetFormulas: universe-self
        """
        els = self.universe.elements.domain() - self.elements.domain()
        dom = self.universe.elements[els]
        if isinstance(self.formula, Not):
            f = self.formula.child
        else:
            f = Not(self)
        name = str(Not(str(self)))
        if self.universe is not None:
            name = self.universe.get_label(els, name)
        return SetFormula(f, dom, self.universe, name)

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

    @Variable.universe.setter
    def universe(self, universe):
        """
        Propagate universe info to subsets

        Args:
            universe (Domain): universe
        """
        self._universe = universe
        if isinstance(self.formula, And) or isinstance(self.formula, Or):
            if isinstance(self.formula.left, SetFormula):
                self.formula.left.universe = universe
            if isinstance(self.formula.right, SetFormula):
                self.formula.right.universe = universe
        elif isinstance(self.formula, Not):
            if isinstance(self.formula.child, SetFormula):
                self.formula.child.universe = universe
        else:
            pass

    def enumerate_elements(self):
        """
        Enumerates the labels corresponding to the multiset

        Returns:
            str: string representation of the multiset
        """
        labels = [f"{self.get_label(e)}" for e in self]
        avg_length = sum([len(lab) for lab in labels]) / len(labels)
        separator = ",\n\t" if avg_length > 7 else ", "
        lab_list = separator.join(labels)
        if avg_length > 7:
            enum = f"{{\n\t{lab_list}\n}}"
        else:
            enum = f"{{{lab_list}}}"
        return enum

    # def take(self, n):
    #     taken = Multiset.take(n)
    #     return SetFormula(
    #         self.formula,
    #         taken.elements,
    #         self.universe,
    #         name=taken.name,
    #         labels=self.labels,
    #     )


class Universe(SetFormula):
    """
    A SetFormula aware of being the Universe. Keeps track of the user-defined labels associated to individual entities and sets
    """

    def __init__(
        self, formula, elems, name="universe", labels_entity={}, labels_set={}
    ):
        super().__init__(formula, elems, None, name)
        self.labels_entity = labels_entity
        self.labels_set = labels_set
        if len(support(elems)) > 0:
            self.labels_set[support(elems)] = self.name

    def add_set_label(self, domain, label):
        self.labels_set[domain] = label
        complement = support(self.elements) - domain
        if complement not in self.labels_set:
            self.labels_set[complement] = str(Not(label))

    def get_label(self, key, default="???"):
        if isinstance(key, int):
            return self.labels_entity.get(key, default)
        elif isinstance(key, P.Interval):
            if key.empty:
                return "empty"
            self.labels_set.setdefault(key, default)
            return self.labels_set[key]
        else:
            raise Exception(f"Unknown key type for label: {type(key)}")

    def get_labels(self):
        return self.labels_entity, self.labels_set

    # def find_label(self, key):
    #     """Approximates a smart way to get short labels for parts of the universe combining known labels
    #         But it is not guaranteed that the lable is the best possible

    #     Args:
    #         key (Interval): the interval(s) corresponding to the elements of the part

    #     Returns:
    #         str: a human-readable label for the given set
    #     """
    # parents = []
    # for domain in self.labels_set.keys():
    #     if domain.contains(key):
    #         parents.append(domain)
    # if len(parents) == 0:
    #     children = []
    #     for domain in self.labels_set.keys():
    #         if key.contains(domain):
    #             children.append(domain)
    # min_parent = P.empty()
    # for p in parents:
    #     if len(p) < len(min_parent):
    #         min_parent = p
    # print(min_parent, key)
    # for domain in self.labels_set.keys():
    #     if min_parent - domain == key:
    #         lab_mp = self.labels_set[min_parent]
    #         lab_dom = self.labels_set[domain]
    #         lab_dom = lab_dom[1:] if lab_dom.startswith("¬") else f"¬{lab_dom}"
    #         lab_key = f"{lab_mp} ∧ {lab_dom}"
    #         self.labels_entity
    #         return lab_key
    # enum = [self.get_label(k) for k in P.iterate(key, step=1)]
    # lab_list = ", ".join(enum)
    # return f"{{{lab_list}}}"

    def copy(self):
        return Universe(self.formula, self.elements, name=self.name)

    @property
    def universe(self):
        return self
