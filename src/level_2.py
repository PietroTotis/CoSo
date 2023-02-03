import portion as P

from .configuration import Variable
from .util import *


class LiftedSet(Variable):
    """
    Represents a lifted subset of the universe/a set

    Attributes
    ----------
    name : str
        not important at the moment
    size : SizeFormula
        describes the set of cardinality values that are valid for this set
    ccs : [CCountings]
        describes properties of the set in terms of counting formulas
    relevant : [CCounting]
        defines which groups of indistinguishable objects need to be accounted for
    source : SetFormula
        the set of which the LiftedSet is subset
    """

    def __init__(self, universe, size, ccs=[]):
        """
        A LiftedSet is a set of constraints restricting the powerset of the universe.

        Args:
            universe (Multiset): universe of the problem
            size (CSize): constrains the size of the set (at least non empty)
            ccs ([CConstraints], optional): further counting constraints. Defaults to [].
        """
        super().__init__(universe, f"part. of {universe}")
        self.universe = universe
        self.size = size
        self.ccs = self.compact_ccs(ccs)
        self.histogram = {}
        self.check_bound()

    def __and__(self, rhs):
        """
        Conjoin with another set of constraints

        Args:
            rhs (LiftedSet): other

        Returns:
            LiftedSet: conjunction of constraints
        """
        # name = f"{self.name} /\ {rhs.name}"
        size = self.size & rhs.size
        ccs = self.ccs + rhs.ccs
        ccs = self.compact_ccs(ccs)
        return LiftedSet(self.universe, size, ccs)

    def __contains__(self, constr):
        """
        Check if a constraint is present in the lifted set

        Args:
            constr (Constraint): Counting or CSize

        Returns:
            bool: true if there is a constraint implying constr
        """
        if hasattr(constr, "formula"):
            for cc in self.ccs:
                if cc.formula == constr.formula:
                    return cc.values in constr.values
            return False
        else:
            return self.size.values in constr.values

    def __eq__(self, rhs):
        """
        Check if two lifted sets have the same constraints

        Args:
            rhs (LiftedSet): other

        Returns:
            bool:
        """
        if self.size != rhs.size:
            return False
        elif self.ccs != rhs.ccs:
            return False
        elif self.histogram != rhs.histogram:
            return False
        else:
            return True

    __hash__ = Variable.__hash__

    def __str__(self):
        s = f"{self.name}: {self.size}"
        if len(self.ccs) > 0:
            s += "\n"
            s += "\t counting: "
            s += ",".join([str(c) for c in self.ccs])
        if len(self.histogram) > 0:
            s += "\n\t histogram:" + str(self.histogram)
        return s

    def add_cc(self, cc):
        """
        Add a new counting constraint replacing infs with int numbers

        Args:
            cc (CConstraint): new counting constraint

        Returns:
            LiftedSet: return a new lifted set wit the additional constraint
        """
        s = self.size.values.upper
        s = s - 1 if self.size.values.right == P.OPEN else s
        if cc.values.upper == P.inf:
            safe_vals = cc.values.replace(upper=s + 1)
        else:
            safe_vals = cc.values
        if cc.values.lower == P.inf:
            safe_vals = safe_vals.replace(lower=s)
        cc.values = safe_vals
        ccs = self.ccs + [cc]
        return LiftedSet(self.universe, self.size, ccs)

    def bound(self, rv_set):
        """
        Check if the number of entities from a relevant set is determined by constraints

        Args:
            rv_set (Multiset): relevant set

        Returns:
            bool: wether the number of entities from rv_set can be inferred from the constraints
        """
        bound = False
        if rv_set == self.universe:
            bound = is_singleton(self.size)
        else:
            for cc in self.ccs:
                if cc.formula == rv_set and is_singleton(cc.values):
                    self.histogram[rv_set] = cc.values.lower
                else:
                    bound = bound or is_singleton(cc.values)
        return bound

    def compact_ccs(self, counts):
        """
        Clean the list of counting constraints removing trivial/duplicates

        Args:
            counts ([CConstraints]): list of counting constraints

        Returns:
            [CConstraints]: compacted list
        """
        compact = []
        remove = []
        indexes = range(0, len(counts))
        for i in indexes:
            cc1 = counts[i]
            for j in [j for j in indexes if j > i]:
                cc2 = counts[j]
                if cc1 != cc2:
                    if cc1.formula == cc2.formula:
                        new_interval = cc1.values & cc2.values
                        merged_cc = cc1.copy()
                        merged_cc.values = new_interval
                        compact.append(merged_cc)
                        remove += [i, j]
                else:  # cc1 == cc2:
                    remove += [i]
        keep = [i for i in indexes if not i in remove]
        compact += [counts[i] for i in keep]
        final = len(remove) == 0
        result = compact if final else self.compact_ccs(compact)
        return result

    def check_bound(self):
        """
        Removes inf upper bound from counting constraints
        """
        for cc in self.ccs:
            if cc.values.upper == P.inf:
                ub = self.size.values.upper + 1
                max_int = P.closed(0, ub)
                cc.values = cc.values & max_int

    def copy(self):
        size = self.size.copy()
        ccs = self.ccs.copy()
        hist = self.histogram.copy()
        ls = LiftedSet(self.universe, size, ccs)
        ls.histogram = hist
        return ls

    def disjoint(self, constr):
        """
        Check if a constraint is incompatible with current constraints

        Args:
            constr (Constraint): CConstraint or CSize

        Returns:
            bool: test result
        """
        if hasattr(constr, "formula"):  # counting formula
            dis_ccs = False
            for cc in self.ccs:
                if cc.formula == constr.formula:
                    if cc.values & constr.values == P.empty():
                        dis_ccs = True
            return dis_ccs
        else:  # size formula
            disj = (self.size.values & constr.values) == P.empty()
            return disj

    def feasible(self, rv_set, n, n_class=1):
        """Check if there is some incompatible constraint with the distribution of n entities of
        rv_set over n_class exchangeable variables

        Args:
            rv_set (SetFormula): the domain formula being distributed
            n (int): number of entities
            n_class (int, optional): number of exchangeable variables together with self. Defaults to 1.

        Returns:
            Boolean: True if no constraint is violated by the distribution
        """
        # print("-----")
        # print(self, rv_set, n, n_class)
        if rv_set == self.universe:
            # we are distributing the universe: a size constraint
            size_class = self.size.values.replace(
                upper=lambda v: n_class * v, lower=lambda v: n_class * v
            )
            feasible = n in size_class
            # if not feasible:
            # print(f"Unfeasible size: {n} not in {size_class}")
            return feasible
        else:
            for cc in self.ccs:
                if rv_set in cc.formula:
                    cc_class = cc.values.replace(
                        upper=lambda v: n_class * v, lower=lambda v: n_class * v
                    )
                    # we would have too many elements for a counting constraint
                    if cc_class.empty:
                        return False
                    lb = (
                        cc_class.lower
                        if cc_class.left == P.CLOSED
                        else cc_class.lower + 1
                    )
                    ub = (
                        cc_class.upper
                        if cc_class.right == P.CLOSED
                        else cc_class.upper - 1
                    )
                    feasible = n <= ub
                    if not feasible:
                        # print(f"Unfeasible because {n} {rv_set} > {ub} {cc.formula}")
                        return feasible
                    # now check histograms: because rv_set in cc.formula we start with
                    # n entities that satisfy cc.formula
                    n_cc = 0
                    unfixed = []
                    # then we count which entities are fixed that satisfy cc.formula
                    for rvs in self.histogram:
                        if rvs in cc.formula and rvs != rv_set:
                            if self.histogram[rvs] == -1:
                                unfixed.append(rvs)
                            else:
                                n_cc += self.histogram[rvs]
                    if len(unfixed) == 0:
                        # if everything is fixed then the constraint should be satisfied
                        if n_cc + n not in cc_class:
                            # print(f"Everything related to {cc} is fixed but still unsat")
                            return False
                    else:
                        # if the sizes of the unfixed (in cc.formula) are too small to satisfy
                        # a lower bound of the constraint then unsat
                        sizes = sum([rvs.size() for rvs in unfixed])
                        if (sizes + n_cc + n) * n_class < lb:
                            # print(f"{sizes} entities from other relevant sets and {n_class} variables are not enough for {cc}")
                            return False
            if len(self.histogram) > 1:
                # check lasts unfixed value for histogram: everything not yet sat should become sat
                vals = [n for n in self.histogram.values() if n > -1]
                if (
                    len(vals) == len(self.histogram) - 1
                    and self.histogram[rv_set] == -1
                ):  # last relevant set: check size
                    fixed = sum(vals)  # respect overall size
                    if fixed + (n / n_class) not in self.size.values:
                        # print(f"{rv_set} is the last to be fixed but {fixed+(n/n_class)} is not a valid size")
                        return False
            return True

    def set_labels(self, labels):
        super().set_labels(labels)
        for cc in self.ccs:
            cc.set_labels(labels)

    def rv_size(self, relevant):
        """
        Return the (possible) number of entities from a specific relevant subset

        Args:
            relevant (SetFormula): a relevant set

        Returns:
            portion: interval of allowed size of the subset of entities in relevant
        """
        for cc in self.ccs:
            if cc.formula == relevant:
                return cc.values

    def size_is_defined(self):
        """
        Check if the possible sizes are a unique value

        Returns:
            bool: result
        """
        s = 0
        for e in self.size.values:
            if not e.empty:
                if e.left == P.CLOSED and e.right == P.CLOSED:
                    s += e.upper - e.lower + 1
                elif e.left == P.OPEN and e.right == P.OPEN:
                    s += e.upper - e.lower - 1
                else:
                    s += e.upper - e.lower
        return s == 1

    def relevant(self):
        """
        Returns the sets that are relevant for this LiftedSet

        Returns:
            {Multiset}: the relevant sets
        """
        relevant = set([cc.formula for cc in self.ccs])
        if len(relevant) == 0:
            relevant = {self.universe}
        return relevant

    def satisfies(self, constraint, lb, ub):
        """
        Check if a constraint is already implied by others

        Args:
            constraint (CConstraint): constraint to be checked
            lb (int): lower bound on how many sets have to satisfy the level 1 property
            ub (int): upper bound on how many sets have to satisfy the level 1 property

        Returns:
            bool: result
        """
        sat = None
        for cc in self.ccs:
            if cc.formula in constraint.formula:
                int = constraint.values & cc.values
                if cc.values in constraint.values:
                    sat = True
                elif int.empty:
                    sat = False
        s = constraint.formula.size()
        cc_lb, cc_ub = interval_closed(
            constraint.values, ub_default=self.size.values.upper
        )
        if cc_lb > cc_ub:
            # asking to have more elements than the size of the part
            sat = False
        if cc_lb * lb > s:
            # asking to distribute more elements than available
            sat = False
        return sat
