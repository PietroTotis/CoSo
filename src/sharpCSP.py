import portion as Int
import math
import itertools

from .configuration import *
from .count import *
from .level_1 import *
from .level_2 import *
from .logger import ProblemLog


class Unsatisfiable(Exception):
    def __init__(self, value, log=None):
        self.value = value
        self.log = log
        if self.log is not None:
            self.log.solution = Zero()

    def __str__(self):
        return repr(self.value)


class SharpCSP(object):
    """
    Counts the possible assignments to a set of vars given choice/count constraints

    Attributes:
    vars : [Variable]
        the variables representing the problem configuration
    config : Configuration
        configuration properties associated to vars
    # agg_f : [AggFormula]
    #     aggrgate constraints on integers
    choice_f : [CPosition]
        choice formulas to enforce
    count_f : [CCounting]
        set of size constraints to enforce
    fixed_choices: int
        subsets only: number of variables already set by choice constraints
    lvl: int
        nesting level for logging
    """

    def __init__(
        self,
        vars,
        config,
        choice_f,
        count_f,
        universe=None,
        lvl=0,
        caption="Subproblem",
        id="1",
        debug=True,
    ):
        # self.agg_f = agg_f
        self.choice_f = choice_f
        self.count_f = count_f
        self.debug = debug
        self.fixed_choices = 0
        self.log = ProblemLog(
            vars,
            config,
            choice_f,
            count_f,
            universe=universe,
            level=lvl,
            id=id,
            caption=caption,
            debug=debug,
        )
        self.id = id
        self.lvl = lvl
        self.n_vars = len(vars)
        self.config = config
        self.vars = vars
        self.universe = universe

    def __str__(self):
        str = "----------\n"
        str += f"Solving {self.n_vars} vars:\n"
        for v in self.vars:
            str += f"\t{v}\n"
        str += "Choice constraints:\n"
        for c in self.choice_f:
            str += f"\t{c}\n"
        str += "Count constraints:\n"
        for c in self.count_f:
            str += f"\t{c}\n"
        # str += "Aggregate constraints:\n"
        # for c in self.agg_f:
        #     str += f"\t{c}\n"
        str += "----------"
        return str

    # def apply_aggs(self):
    #     """
    #     either propagates or counts an aggregate constraint on integers
    #     """
    #     agg = self.agg_f[0]
    #     others = self.agg_f[1:]
    #     if agg.op == min:
    #         return self.propagate_min(agg, others)
    #     else:
    #         return self.propagate_max(agg, others)

    def apply_choice(self, chf):
        """
        Sets variables according to the choice constraints
        """
        self.log.add_relevant_set(chf.formula)
        if isinstance(chf, CPosition):
            if chf.pos - 1 < self.n_vars:
                self.vars[chf.pos - 1] = self.vars[chf.pos - 1] & chf.formula
                # self.debug("Choice set: ", self.vars[chf.pos - 1])
                self.log.propagation(chf)
            else:
                # self.debug(f"Ignoring {chf} because out of range ({self.n_vars})")
                self.log.warning(f"Ignoring {chf} because out of range ({self.n_vars})")

    def apply_count_constraint(self, cc, others):
        """
        Applies the constraints to satisfy 'current' count: for each admissible value n (lb<=n<=ub)
        of the property p in cc, constrain the problem to n variables satisfying p and sum over the
        different (independent) counts

        Args:
            cc (CCounting): counting formula to apply
            others ([CCounting]): remaining counting formulas

        Returns:
            Count: count of the propagation
        """
        # self.debug("Propagating ", cc)
        al = self.log.action(f"Considering constraint {cc}")
        self.log.add_relevant_set(cc.formula)
        out_values = Int.closedopen(0, self.config.size.values.upper) - cc.values
        if cc.values.atomic:
            return self.apply_simple_interval_count(cc, others, al)
        elif out_values.atomic:
            # ignore and then remove the unsat cases
            # self.debug(f"Relax {cc} and remove unsat")
            self.log.detail(al, f"Relax {cc} and remove unsat (={out_values})")
            count_ignore = self.split_on_constraints(1, others, op="add")
            not_cc = CCounting(cc.formula, out_values)
            count_not = self.split_on_constraints(
                1, [not_cc] + others, op="sub", id="2"
            )
            count = count_ignore - count_not
        else:
            # there are multiple separate intervals
            self.log.detail(al, f"Breaking intervals into equalities")
            count = Zero()
            for interval in cc.values:
                partial_cc = CCounting(cc.formula, interval)
                partial_count = self.solve_subproblem(
                    self.vars,
                    self.config,
                    {},
                    [partial_cc] + others,
                    caption=f"Case with {interval} {cc.formula}",
                    op="add",
                    id=str(interval),
                )
                # partial_count = self.apply_count_constraint(partial_cc, others)
                # self.log.add_subproblem("add", partial_count.log)
                count += partial_count
            return count

    def apply_constraints(self):
        """
        If vars are exchangeable we can apply a count otherwise we need to split and consider the combinations of count constraints

        Returns:
            Count: problem sollution
        """
        ex_classes = self.exchangeable_classes()
        if self.config.lvl1():
            if len(ex_classes) > 1:
                count = self.count_non_exchangeable(ex_classes)
            else:
                count = self.count_exchangeable()
        else:
            count = self.shatter_partitions()
        return count

    def apply_inverse_count(self, cc, others, out_values, action):
        """
        If there are more value included than excluded
        ignore and then remove the unsat cases

        Args:
            cc (CCounting): the constraint to be propagated
            others ([Ccounting]): rest of constraints to propagate later
            out_values (Interval): values excluded by the constraint
            action (ActionLog): logging

        Returns:
            Count: nr. solutions with this constraint propagated
        """
        # self.debug(f"Relax {cc} and remove unsat")
        self.log.detail(action, f"Relax {cc} and remove unsat")
        count_ignore = self.split_on_constraints(1, others, op="add")
        not_cc = CCounting(cc.formula, out_values)
        count_not = self.split_on_constraints(1, [not_cc] + others, op="sub", id="2")
        return count_ignore - count_not

    def apply_propositional_count(self, cc, others, ub, action):
        """
        Break interval into equalities with propositionalization of each i in the interval

        Args:
            cc (CCounting): the constraint to be propagated
            others ([Ccounting]): rest of constraints to propagate later
            ub (int): Upper bound from configuration
            action (ActionLog): logging

        Returns:
            Count: nr. solutions with this constraint propagated
        """
        # self.debug(f"Expanding bounds {cc.values}...")
        values = cc.values.replace(upper=ub, right=Int.CLOSED)
        self.log.detail(
            action,
            f"Shatter values {cc.values}: {[i for i in Int.iterate(values, step=1)]}",
        )
        count = Zero()
        # optimization:if we expand values and ...=n is not satisfiable
        # then we can stop beacause ...>= n must be unsat from n on
        stop_at_unsat = False
        for i in Int.iterate(values, step=1):
            if not stop_at_unsat:
                cc_eq = CCounting(cc.formula, Int.singleton(i))
                # self.debug(f"Propagate lower bound #{cc.formula}={lb}...")
                # self.log.detail(al, f"Propagate lower bound #{cc.formula}={lb}...")
                count_case = self.solve_subproblem(
                    self.vars,
                    self.config,
                    {},
                    [cc_eq] + others,
                    caption=f"Case with {i} {cc.formula}",
                    op="add",
                    id=str(i),
                )
                if count_case.is_zero():
                    self.log.detail(action, f"Stopping at {i} which is already 0")
                    stop_at_unsat = True
                count += count_case
        return count

    def apply_simple_count(self, cc, others, action):
        """
        What everything boils down to: propagate the constraint that exactly n
        objects/variables have a property

        Args:
            cc (CCounting): the constraint to be propagated
            others ([Ccounting]): rest of constraints to propagate later
            action (ActionLog): logging stuff

        Returns:
            Count: nr. solutions with this constraint propagated
        """
        try:
            # self.log.detail(al, f"Propagating {cc}")
            choices, prop_vars = self.propagate_cc(cc)
            self.vars = prop_vars
            count = self.split_on_constraints(choices, others)
        except Unsatisfiable as u:
            # self.log.add_subproblem("unsat", self.log)
            # print(self.log)
            count = Zero(tip=u.value)
        return count

    def apply_simple_interval_count(self, cc, others, action):
        """This function collects a few optimizations we can do to count on
            a continuous interval.
        Args:
            cc (CCounting): the constraint with continuous values
            others ([CCounting]): rest of constraints to be propagate
            action (ActionLog): logging stuff

        Returns:
            Count : the count corresponding to propagating this constraint
        """
        max_entities = min(self.n_vars, cc.formula.size())
        # atomic means that there is a single continuous interval
        lb, ub = interval_closed(cc.values, ub_default=max_entities)
        all_values = Int.closed(0, self.n_vars)
        out_values = all_values - cc.values
        if lb == ub:
            return self.apply_simple_count(cc, others, action)
        elif self.config.type != "subset" and (ub >= self.n_vars):
            self.log.detail(
                action,
                f"Fix at least ({cc.values.left}) {cc.values.lower} {cc.formula} and let rest unconstrained",
            )
            return self.apply_simple_count(cc, others, action)
        elif lb == 0 and ub >= self.n_vars:
            # trivial
            return self.split_on_constraints(1, others)
        elif lb == 0 and ub < self.n_vars:
            # propagate n_vars-ub ¬cc.formula and leave rest unconstrained
            ub = self.n_vars - ub
            # self.debug(f"Invert {cc}: at least {ub} {cc.formula.neg()}")
            self.log.detail(action, f"Invert {cc}: at least {ub} {cc.formula.neg()}")
            interval_not = Int.closedopen(ub, Int.inf)
            cc_not = CCounting(cc.formula.neg(), interval_not)
            return self.apply_count_constraint(cc_not, others)
        elif len(out_values) < len(cc.values):
            return self.apply_inverse_count(cc, others, out_values, action)
        else:
            return self.apply_propositional_count(cc, others, ub, action)

    def choose_cc(self):
        """
        Heuristics to decide which counting constraint to propagate when variables are exchangeable

        Returns:
            (CCounting, [CCounting]): the chosen counting constraint to be propagated and the list of the others
        """
        choice = self.count_f[0]
        others = self.count_f[1:]
        ex_class = self.vars[0]
        if isinstance(ex_class, SetFormula):
            for cc in self.count_f:
                same_class = cc.formula == ex_class
                max = cc.values.lower == self.n_vars
                zero = cc.values.upper == 0
                if same_class or max or zero:
                    choice = cc
                    others = self.count_f
                    others.remove(cc)
                    return choice, others
            return choice, others
        else:
            return choice, others

    def compact_ccs(self, counts):
        """
        Simplify some counting constraints, i.e. turn (#p > n), (#p < m) into (#p in [n,m])

        Args:
            counts ([CCounting]): list of counting constraints to be simplified

        Raises:
            Unsatisfiable: if we find some incompatible constraints then counts is unsatisfiable

        Returns:
            [CCounting]: equivalent list of counting constraints compacted
        """
        compact = []
        indexes = range(0, len(counts))
        remove = []
        for i in indexes:
            cc1 = counts[i]
            if cc1.values.empty:
                raise Unsatisfiable("Empty set of valid counts", self.log)
            satisfied, not_satisfied, maybe = self.count_satisfied(cc1)
            if (
                len(maybe) == 0 and len(satisfied) in cc1.values
            ):  # remove ccs already satisfied
                remove += [i]
            else:
                for j in [j for j in indexes if j > i]:
                    cc2 = counts[j]
                    if cc1 != cc2:
                        if cc1.formula == cc2.formula:  # combination is not satisfiable
                            new_int = cc1.values & cc2.values
                            if new_int.empty:
                                raise Unsatisfiable("Conflicting constraints", self.log)
                            else:  # combination satisfiable: merge intervals into single cc
                                merged_cc = CCounting(cc1.formula, new_int)
                                compact.append(merged_cc)
                                remove += [i, j]
                    else:  # cc1 == cc2:
                        # opposite formulas both 0
                        remove += [i]
        keep = [i for i in indexes if not i in remove]
        compact += list(map(lambda i: counts[i], keep))
        final = len(remove) == 0
        result = compact if final else self.compact_ccs(compact)
        return result

    def count_exchangeable(self):
        """
        Take one count constraint and apply it

        Returns:
            Count: the count of solutions of the problem
        """
        # self.debug("Counting exchangeable...")
        cc, others = self.choose_cc()
        count = self.apply_count_constraint(cc, others)
        return count

    def count_non_exchangeable(self, ex_classes):
        """
        If vars are not exchangeable we split the problem into class | rest
        Counting constraints need to be split according to the possible combinations
        Example:
        Split classes french | students with ccs on dutch=2 and french=1
        for each cc prop==n the count can be split between the two classes with i elements on the left and n-i elements on the right, i.e. french==1 | french==0 and french==0 | french==1.
        But with many ccs we need to consider all the combinations of the two so:
        french==0 & dutch==0 | french==1 & dutch==2
        french==1 & dutch==0 | french==0 & dutch==2 ...

        Args:
            ex_classes ([SetFormulas/LiftedSet]): variables grouped by exchangeable classes

        Returns:
            Count: count of the solutions of the problem
        """
        # TODO: optimize by propagating first constraints on sets that are disjoint from others
        split_class_vars, rest_classes_vars = self.split_ex_classes(ex_classes)
        combs_split_class, combs_rest_classes = self.shatter_count_formulas(
            split_class_vars, rest_classes_vars
        )
        al = self.log.action("Shattering constraints")
        tot_count = Zero()
        for i in range(0, len(combs_rest_classes)):
            try:
                comb_split_class = self.compact_ccs(combs_split_class[i])
                comb_rest_classes = self.compact_ccs(combs_rest_classes[i])
                # self.debug(
                #     f"Solving combination {i}: {comb_split_class} // {comb_rest_classes}"
                # )
                self.log.shatter_case(al, i, comb_split_class, comb_rest_classes)
                split_args = [
                    split_class_vars,
                    rest_classes_vars,
                    list(comb_split_class),
                    list(comb_rest_classes),
                ]
                if self.config.lvl1() and not self.config.with_repetition():
                    count = self.split_inj(*split_args)
                elif self.config.lvl1() and self.config.with_repetition():
                    count = self.split(*split_args, shatter_id=i)
                else:
                    count = self.split_partitions(*split_args)
                desc_csc = ", ".join([str(c) for c in comb_split_class])
                desc_crc = ", ".join([str(c) for c in comb_rest_classes])
                self.log.detail(al, f"Split {i+1}: {desc_csc} // {desc_crc}")
                # self.debug(f"Split combination count: {count}")
            except Unsatisfiable as u:
                count = Zero(tip=u.value)
            tot_count += count
        # self.debug(f"Shatter count: {tot_count}")
        return tot_count

    def shatter_count_formulas(self, split_class_vars, rest_classes_vars, ccs=None):
        """
        Returns a shattering of the counting constraints for the given split.
        Given a split, for each size constraints compute the set of splits that shatter the constraint.
        Then combine each different cc shattering with the other constraints' possible splits (cartesian product).
        """
        # self.debug("Shattering size constraints...")
        # self.log.action("Shattering")
        count_f = self.count_f if ccs is None else ccs
        n_split = len(split_class_vars)
        n_rest = self.n_vars - n_split
        ccs_split_class = []
        ccs_rest_classes = []
        for cc in count_f:
            # cases_split_class = []  # split for cc (left)
            # cases_rest_classes = [] # split for cc (right)
            # vals = cc.values
            # if cc.values.upper >= self.n_vars:
            #     interval_any = cc.values.replace(upper=n_split, right=Int.CLOSED)    # adjust upper bound to number of vars
            #     any = CCounting(cc.formula, interval_any)    # add case >= n in split class then work on each =i <n
            #     cases_split_class.append(any)
            #     cc_rest_classes = CCounting(cc.formula, Int.closed(0,n_rest))
            #     cases_rest_classes.append(cc_rest_classes)
            #     if cc.values.left == Int.OPEN:
            #         n_cases = cc.values.lower
            #     else:
            #         n_cases = cc.values.lower-1
            # else:
            #     n_cases = cc.values.upper
            # for i in range(0,n_cases+1):
            # for partition in self.get_feasible_splits(split_class_vars[0],n_split, n_rest, cc.formula):
            #     l_value, r_value = partition
            #     cc_split_class = CCounting(cc.formula, Int.singleton(l_value))
            #     cc_rest_classes = CCounting(cc.formula, Int.singleton(r_value))
            #     # if self.is_feasible_split(split_class_vars[0],n_split,n_rest,cc_split_class,cc_rest_classes):
            #     cases_split_class.append(cc_split_class)
            #     cases_rest_classes.append(cc_rest_classes)
            cases_split_class, cases_rest_classes = self.get_feasible_splits(
                split_class_vars[0], n_split, n_rest, cc
            )
            ccs_split_class.append(cases_split_class)
            ccs_rest_classes.append(cases_rest_classes)
        combs_split_class = list(itertools.product(*ccs_split_class))
        combs_rest_classes = list(itertools.product(*ccs_rest_classes))
        return combs_split_class, combs_rest_classes

    def count(self, vars=None):
        if self.config.lvl1() and self.config.labelled():
            count = self.count_sequence(vars)
        elif self.config.lvl1() and not self.config.labelled():
            count = self.count_subsets(vars)
        elif self.config.lvl2():
            count = self.count_partitions(vars)
        else:
            raise Exception("Unknown configuration type")
        return count

    def count_satisfied(self, cc, var_list=None):
        """Checks for each variable if the property counted

        Args:
            property (SetFormula): _description_
            var_list ([SetFormula]/[LiftedSet], optional): _description_. Defaults to None.
            n (int, optional): how many exchangeable variables need to have the property,
                               useful for checking second level properties. Defaults to 1.

        Returns:
            ([int],[int],[int]): sat/maybe/unsat lists of indexes of variables in each group
        """
        vars = var_list if var_list is not None else self.vars
        property = cc.formula
        sat = []
        not_sat = []
        maybe = []
        if isinstance(property, SetFormula):
            for i, v in enumerate(vars):
                if v in property:
                    sat.append(i)
                elif v.disjoint(property):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        elif isinstance(property, CSize):
            for i, v in enumerate(vars):
                if v.size in property:
                    sat.append(i)
                elif v.size & property == CSize("", Int.empty()):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        else:
            lb, ub = interval_closed(cc.values, ub_default=len(vars))
            for i, v in enumerate(vars):
                satisfies = v.satisfies(property, lb, ub)
                if satisfies is None:
                    maybe.append(i)
                elif satisfies:
                    sat.append(i)
                else:
                    not_sat.append(i)
        return (sat, not_sat, maybe)

    def disjoint(self, classes):
        disj = True
        vars = list(classes.keys())
        indexes = range(0, len(vars))
        for c1 in indexes:
            for c2 in [j for j in indexes if j > c1]:
                disj = disj and vars[c1].disjoint(vars[c2])
        return disj

    def exchangeable_classes(self, var_list=None):
        """
        Groups exchangeable variables in a dict

        Args:
            var_list ([SetFormula/LiftedSet], optional): _description_. Defaults to None.

        Returns:
            {}: _description_
            TODO
        """
        vars = var_list if var_list is not None else self.vars
        classes = {}
        for v in vars:
            new = True
            for c in classes:
                if v == c:
                    classes[c].append(v)
                    new = False
            if new:
                classes[v] = [v]
        return classes

    def filter_domain(self, cases, dformula):
        """
        Given a domain formula dformula and the domain formula case_f of n elements shrinks dformula by (up to) n elements satisfying case_f
        """
        for case in cases:
            n = cases[case]
            inter = dformula & case
            exclude = inter.take(n)
            if exclude.size() > 0:
                df = dformula - exclude
        return df

    def filter_domains(self, cases, variables):
        # self.debug("Filtering domains...")
        al = self.log.action("Filtering domains")
        dfs = variables
        for case in cases:
            # self.debug(f"  Case {case}")
            self.log.detail(al, f"Case {case}")
            n = cases[case]
            if n > 0:
                pool = case
                for df in dfs[1:]:
                    if not df.disjoint(case):
                        pool = pool & df
                exclude = pool.take(n)
                if exclude.size() > 0:
                    updated_dfs = []
                    for df in dfs:
                        # self.debug(f"  Filtering {n} {case} out of {df}:")
                        self.log.detail(al, f"Filtering {n} {case} out of {df}")
                        df -= exclude
                        updated_dfs.append(df)
                    dfs = updated_dfs
                else:
                    dfs = [exclude] * n
        return dfs

    def get_vars(self, indexes):
        return [self.vars[v] for v in indexes]

    def get_subproblem_universe(self, vars, given):
        """Return the correct universe depending on the given one

        Args:
            var (Variable): a variable of the configuration
            given (MultiSet): optional
        """
        ex_classes = self.exchangeable_classes(var_list=vars)
        if given is not None:
            # Promote the multiset to Universe
            u = Universe(
                given.formula,
                given.elements,
                given.name,
                *self.universe.get_labels(),
            )
        elif len(ex_classes) == 1:
            u = Universe(
                vars[0].formula,
                vars[0].elements,
                vars[0].name,
                *self.universe.get_labels(),
            )
        else:
            u = self.universe
        return u

    def histogram(self, var_list=None):
        """
        Returns each exchangeable class with the number of variables belonging to the class
        """
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        for c in ex_classes:
            ex_classes[c] = len(ex_classes[c])
        return ex_classes

    def integer_partitions(self, n):
        """
        From http://jeromekelleher.net/generating-integer-partitions.html
        """
        a = [0 for i in range(n + 1)]
        k = 1
        y = n - 1
        while k != 0:
            x = a[k - 1] + 1
            k -= 1
            while 2 * x <= y:
                a[k] = x
                y -= x
                k += 1
            l = k + 1
            while x <= y:
                a[k] = x
                a[l] = y
                yield a[: k + 2]
                x += 1
                y -= 1
            a[k] = x + y
            y = x + y - 1
            yield a[: k + 1]

    def integer_k_partitions(self, n, k):
        """
        Partitions n in k (possibly 0) integers
        """
        parts = []
        for p in self.integer_partitions(n):
            if len(p) < k:
                padded = p + [0] * (k - len(p))
                parts.append(padded)
            elif len(p) == k:
                parts.append(p)
        return parts

    def get_feasible_splits(self, split_class_var, n_split, n_rest, cc):
        disjoint = split_class_var.disjoint(cc.formula)
        # if not isinstance(cc.formula, CCounting):
        if isinstance(cc.formula, SetFormula):
            included = split_class_var in cc.formula
        else:
            # TODO: do the same check for level 2 configurations
            included = False
        cases_split_class = []  # split for cc (left)
        cases_rest_classes = []  # split for cc (right)
        for inter in cc.values:
            lb, ub = interval_closed(inter, ub_default=n_split + n_rest)
            if included and n_split > ub:
                # unsat
                pass
            elif disjoint:
                # rest does all
                if n_rest >= lb:
                    cc_split_class = CCounting(cc.formula, Int.singleton(0))
                    cc_rest_classes = CCounting(cc.formula, Int.closed(lb, n_rest))
                    cases_split_class.append(cc_split_class)
                    cases_rest_classes.append(cc_rest_classes)
                else:
                    # unsat
                    pass
            else:
                if included:
                    lb = max(lb, n_split)  # consider only cases where left>=n_split
                for i in range(lb, ub + 1):
                    partitions = self.integer_k_partitions(i, 2)
                    distributions = partitions + [
                        [r, l] for l, r in partitions if l != r
                    ]
                    for l_val, r_val in distributions:
                        if l_val <= n_split and r_val <= n_rest:
                            cc_split_class = CCounting(cc.formula, Int.singleton(l_val))
                            cc_rest_classes = CCounting(
                                cc.formula, Int.singleton(r_val)
                            )
                            cases_split_class.append(cc_split_class)
                            cases_rest_classes.append(cc_rest_classes)
        return cases_split_class, cases_rest_classes

    def is_feasible_split(self, split_class_var, n_split, n_rest, scc, rcc):
        """
        Checks if we ask to observe more properties than available variables
        Or things that split_class_var cannot satisfy
        """
        self.log.action(
            f"Checking feasibility of splitting {split_class_var} as {scc}//{rcc}"
        )
        if scc.values.empty or rcc.values.empty:
            return False
        slb, sub = interval_closed(scc.values, ub_default=n_split)
        rlb, rub = interval_closed(rcc.values, ub_default=n_rest)
        if slb > n_split:
            # self.debug(f"{scc}//{rcc} is unfeasible because not enough vars in split")
            self.log.detail(
                f"{scc}//{rcc} is unfeasible because not enough vars in split"
            )
            return False
        if rlb > n_rest:
            # self.debug(f"{scc}//{rcc} is unfeasible because not enough vars in rest")
            self.log.detail(
                f"{scc}//{rcc} is unfeasible because not enough vars in rest"
            )
            return False
        disjoint = split_class_var.disjoint(scc.formula)
        included = scc.formula in split_class_var
        if disjoint and slb > 0:
            # self.debug(
            #     f"{scc}//{rcc} is unfeasible because split is disjoint from {scc.formula} and {slb} > 0"
            # )
            self.log.detail(
                f"{scc}//{rcc} is unfeasible because split is disjoint from {scc.formula} and {slb} > 0"
            )
            return False
        if included and sub < n_split:
            # self.debug(
            #     f"{scc}//{rcc} is unfeasible because split is included in {scc.formula} and {sub} < {n_split}"
            # )
            self.log.detail(
                f"{scc}//{rcc} is unfeasible because split is included in {scc.formula} and {sub} < {n_split}"
            )
            return False
        return True

    # def debug(self, s, *args):
    #     if self.debug:
    #         strargs = " ".join(list(map(str, args)))
    #         flat = str(s) + strargs
    #         ind = "\t" * self.lvl
    #         lines = flat.split("\n")
    #         indented = list(map(lambda l: ind + l, lines))
    #         final = "\n".join(indented)
    #         print(final)

    def propagate(self, var, property):
        if isinstance(property, SetFormula):
            return var & property
        elif isinstance(property, CSize):
            return var & LiftedSet(self.universe, property)
        elif isinstance(property, CCounting):
            return var.add_cc(property)
        else:
            raise Exception(
                f"unexpected property type {property}: {type(property)}", self.log
            )

    def propagate_cc(self, cc, var_list=None):
        """
        Given a set of exchangeable variables, propagate one size constraint.
        - Check how many vars already satisfy, maybe satisfy, cannot satisfy the property
        - If there is no var that maybe satisfy, either is already sat or not
        - If there are vars that maybe satisfy, set m vars to satisfy the property where m is how many more entities there need to be, if set the others to not satisfy.
        return the number of exchangeable choices the propagation requires
        """
        al = self.log.propagation(cc)
        vars = var_list if var_list is not None else self.vars
        lb, ub = interval_closed(cc.values, ub_default=len(vars))
        atleast = ub == len(vars)
        satisfied, not_satisfied, maybe = self.count_satisfied(cc, vars)
        diff = lb - len(satisfied)
        var_maybe = [vars[v] for v in maybe]
        if diff == 0 and len(maybe) == 0:
            # self.debug(cc, " already satisfied")
            self.log.detail(al, f"is already satisfied")
            return 1, vars
        elif len(maybe) == 0:
            # self.debug(cc, " is unsat here")
            self.log.detail(al, "is unsat with these variables")
            raise Unsatisfiable(
                f"Stopped while trying to propagate {cc} (conflicting with other constraints)"
            )
        else:
            v = var_maybe[0]
            # self.debug(f"{len(maybe)} exchangeable constrainable vars: {v}")
            self.log.detail(al, f"{len(maybe)} exchangeable constrainable vars: {v}")
            if not self.config.labelled():
                n_choices = 1
            else:
                n_choices = math.comb(len(maybe), abs(diff))
            sat_f = self.propagate(v, cc.formula)
            not_sat_f = self.propagate(v, cc.formula.neg())
            for i in maybe:
                if diff < 0 and not atleast:
                    vars[i] = not_sat_f
                    diff += 1
                elif diff > 0:
                    vars[i] = sat_f
                    diff -= 1
                else:  # diff=0
                    if not atleast:
                        vars[i] = not_sat_f
            if diff != 0 and not atleast:
                self.log.detail(al, "is unsat with these variables")
                raise Unsatisfiable(f"Cannot satisfy {cc} :(")
            else:
                self.log.detail(al, self.vars)
                return n_choices, vars

    def relevant_cases_intersection(self, universe, rest_classes):
        """
        Computes recursively the intersections of relevant cases/domains
        """
        first = rest_classes[0]
        if first == universe:
            return self.relevant_cases_intersection(universe, rest_classes[1:])
        combinations = [[first, first.neg()]]
        for dom in rest_classes[1:]:
            # if we have other subsets we can ignore the universe since there is no element
            # in not(universe)
            if not dom == universe:
                comb = [dom, dom.neg()]
                combinations.append(comb)
        combinations = list(itertools.product(*combinations))
        cases = []
        for c in combinations:
            dom_base = c[0]
            if len(c) > 1:
                for dom in c[1:]:
                    dom_base = dom_base & dom
            if dom_base.size() > 0:
                cases.append(dom_base)
        return cases

    def relevant_cases(self, split_df, rest_classes):
        """
        Computes the cases for elements of type split_df w.rt. rest_classes. i.e. if we split on french and we have other classes like dutch and italian,
        we need to consider for alldiff the cases where some of the french are both/neither dutch or italian or one of the two
        """
        # self.debug("Computing case combinations of relevant classes:")
        al = self.log.action("Computing case combinations of relevant classes")
        # if we have one other class consider only positive, since n out of m
        # positives is the same as m-n negatives out of m
        if len(rest_classes) == 0:
            res = []
        elif len(rest_classes) == 1:
            first = rest_classes[0]
            if first != self.universe:
                res = [first, first.neg()]
            else:
                res = [split_df]
        else:
            res = self.relevant_cases_intersection(split_df.universe, rest_classes)
        self.log.detail(al, res)
        for r in res:
            self.log.add_relevant_set(r)
        # self.debug(f"\t{res}")
        return res

    def solve(self):
        for c in self.choice_f:
            self.apply_choice(c)
        try:
            self.count_f = self.compact_ccs(self.count_f)
        except Unsatisfiable as u:
            return Zero(tip=u.value)
        if len(self.count_f) > 0:
            count = self.apply_constraints()
        # elif len(self.agg_f) > 0:
        #     count = self.apply_aggs()
        else:
            count = self.count()
        # self.debug("=========")
        # self.debug("Solution:" + str(count))
        return Solution(count, self.log)

    def split(
        self,
        split_class_vars,
        rest_classes_vars,
        split_class_ccs,
        rest_classes_ccs,
        shatter_id=1,
    ):
        scv = [v.copy() for v in split_class_vars]
        rcv = [v.copy() for v in rest_classes_vars]
        # self.debug(f"Split class :")
        al = self.log.action("Splitting")
        split_class_count = self.solve_subproblem(
            scv,
            self.config,
            [],
            split_class_ccs,
            op="split-left",
            # id=f"{id}",
            shatter_id=shatter_id,
            caption="Split class",
        )
        if split_class_count.is_zero():
            self.log.detail("Left split is 0: ignoring right")
            return split_class_count
        else:
            # self.debug("Rest class: ")
            rest_classes_count = self.solve_subproblem(
                rcv,
                self.config,
                [],
                rest_classes_ccs,
                op="split-right",
                # id=f"{id}",
                shatter_id=shatter_id,
                caption="Rest class",
            )
            count = split_class_count * rest_classes_count
            # self.debug("==========")
            return count

    def solve_subproblem(
        self,
        vars,
        config,
        choice_constr,
        count_constr,
        universe=None,
        op="none",
        caption="Subproblem",
        id="1",
        shatter_id=None,
    ):
        """
        Recursive call to the solver

        Args:
            vars ([Variable]): variables initialization
            config (str): description of the configuration config
            choice_constr ([PConstraint]): positional constraints
            count_constr ([CConstraint]): counting constraints

        Returns:
            Count: solution of the subproblem
        """
        # self.debug(f"\tSubproblem ({config}):")
        vars = [v.copy() for v in vars]

        sub_id = f"{self.id}.{id}"
        if op.startswith("split"):
            if op.endswith("left"):
                sub_id += "L"
            else:
                sub_id += "R"

        subproblem = SharpCSP(
            vars,
            config,
            choice_constr,
            count_constr,
            self.get_subproblem_universe(vars, universe),
            self.lvl + 1,
            caption,
            id=sub_id,
            debug=self.debug,
        )
        try:
            sub_sol = subproblem.solve()

            # logging
            if shatter_id is None:
                self.log.add_subproblem(op, sub_sol.log)
            else:
                if op == "split-left":
                    self.log.add_split_left(sub_sol.log, shatter_id)
                else:
                    self.log.add_split_right(sub_sol.log, shatter_id)
            return sub_sol.count
        except Unsatisfiable as u:
            # self.debug("\t==========\n\tUnsat: 0")
            return Zero(tip=u.value)

    def split_on_constraints(self, n_choices, others, id="1"):
        # self.debug("Splitting on other constraints...")
        al = self.log.action("Splitting on other constraints")
        if len(others) == 0:
            # self.debug("... no other constraints")
            self.log.detail(al, "no other constraints")
            count = self.count()
        else:
            count = self.solve_subproblem(
                self.vars,
                self.config,
                [],
                others,
                op="none",
                caption=self.log.caption,
                id=id,
            )
        count.log = self.log
        res = count.with_choices(n_choices)
        return res

    def split_ex_classes(self, ex_classes):
        """
        Groups the exchangeable classes as left | right by taking one as left and uniting the others
        don't pick universe as split class
        """
        split_class = None
        rest_classes = []
        if self.config.lvl2():
            it_excls = iter(ex_classes)
            split_class = next(it_excls)
            for i in it_excls:
                rest_classes = rest_classes + ex_classes[i]
        else:
            for ex_class in ex_classes:
                if split_class is None and ex_class != self.universe:
                    split_class = ex_class
                else:
                    rest_classes = rest_classes + ex_classes[ex_class]
        # it_excls = iter(ex_classes)
        # split_class = next(it_excls)
        # rest_classes = []
        # for i in it_excls:
        #     rest_classes = rest_classes + ex_classes[i]
        return (ex_classes[split_class], rest_classes)

    #######################
    ## Level 1 methods ##
    #######################

    def count_multisubsets_exchangeable(self, var_list=None):
        """
        Select with replacement from a multiset M: that is equivalent to select with repetition from the support of M:
        find number of entities and apply binomial coefficient

        Args:
            var_list ([SetFormulas]], optional): List of variables representing the mset. Defaults to None.

        Returns:
            Count: count
        """
        vars = var_list if var_list is not None else self.vars
        n = len(vars)
        dom = vars[0]
        indist = dom.elements.find(False)
        indist_sizes = self.get_sizes_indistinguishable(indist)
        indist_size = sum(indist_sizes)
        dist_size = dom.size() - indist_size + len(indist_sizes)
        m = dist_size + n - 1
        # self.debug(f"Multisubsets with all exchangeable: ({dist_size+n-1},{n})={count}")
        al = self.log.action(f"Counting multisubsets with all exchangeable")
        self.log.detail(al, f"With binomial coefficient ({m} {n})")
        sol = Binomial(n, m, f"Choose {n} objects out of {m} distinguishable objects")
        return sol

    def count_multisubsets_non_exchangeable(self, var_list=None):
        """If there are non-exchangeable variables in order to exploit exchangeability classes have to be disjoint.
        Otherwise multiplication is introducing some order. So if dijoint apply multiplication otherwise we need to
        refine the variables that are not disjoint with the relevant parts.

        Args:
            var_list ([SetFormula], optional): List of variables representing the mset. Defaults to None.

        Returns:
            Count: count
        """
        # self.debug("Multisubsets with non-exchangeable...")
        al = self.log.action("Multisubsets with non-exchangeable...")
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        if self.disjoint(ex_classes):
            # self.debug("... but disjoint")
            self.log.detail(al, "... but disjoint")
            count = One()
            for exc in ex_classes:
                s = self.count_multisubsets_exchangeable(ex_classes[exc])
                count *= s
                self.log.detail(al, "Product of multisets of each class")
        else:
            # if classes are not disjoint no multiplication rule
            # self.debug("... and not disjoint: disjoining with relevant parts")
            self.log.detail(al, "... and not disjoint: disjoining with relevant parts")
            relevant = self.relevant_cases_intersection(
                self.universe, list(ex_classes.keys())
            )
            unsplit = [ec for ec in ex_classes if ec not in relevant]
            decompose_class = unsplit[0]
            count = self.shannon_indep_multisubsets(decompose_class, relevant, 0)
        return count

    def shannon_indep_multisubsets(
        self, decompose_class, relevant, previous, var_list=None
    ):
        """If we have two non-disjoit classes for example universe U and a domain D and variables D U U U.
        Then multiplication cannot be applied and we need to refine U. So for each variable that is not a relevant
        part (D or ¬D) we specialize it with the options, e.g. D D U U and D ¬D U U. We do it recursively but
        we are introducing an order, e.g. we can get D D D ¬D or D ¬D D D or D D ¬D D and many other exchangeable combinations.
        With previous we keep track of the last id of relevant part used and we make sure we use only later ones for new
        specializations.

        Args:
            decompose_class (SetFormula): the domain/exchangeability class to be specialized
            relevant ([DomainFOrmula]): the list of partitions that we use to specialize the vars in dijoint domains
            previous (int): index of the last
            var_list ([SetFormula], optional): List of variables representing the mset. Defaults to None.

        Returns:
            Count: count
        """
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        self.log.action("Shattering multiset choices")
        count = Zero()
        for i, r in enumerate(
            relevant
        ):  # specialize variable with each different relevant set
            inter = decompose_class & r  #
            n = len(ex_classes[decompose_class])
            if (
                inter.size() > 0 and i >= previous
            ):  # if the specialized variable is empty or we are not following the order skip
                vars = [
                    v for k in ex_classes for v in ex_classes[k] if k != decompose_class
                ]
                others = [decompose_class] * (
                    n - 1
                )  # specialize one variable with this r then the others recursively if necessary
                vars = vars + [inter] + others
                new_ex_classes = self.exchangeable_classes(vars)
                if len(new_ex_classes) == 1:
                    count += self.count_multisubsets_exchangeable(vars)
                elif self.disjoint(new_ex_classes):
                    count += self.count_multisubsets_non_exchangeable(vars)
                else:
                    unsplit = [ec for ec in new_ex_classes if ec not in relevant]
                    next_decompose_class = unsplit[0]
                    count += self.shannon_indep_multisubsets(
                        next_decompose_class, relevant, i, var_list=vars
                    )
        return count

    def count_sequence(self, var_list=None):
        """
        Look at a sequence problem and figure out what to do:
        - Sequence: cartesian product
        - Sequence with alldiff (permutation):
            - Exchangeable: counting rule (factorial)
            - Non-exchangeable: split/shatter
        """
        vars = var_list if var_list is not None else self.vars
        if self.config.type == "sequence":
            # self.debug("Counting sequences:")
            sol = self.count_sequence_any(vars)
        else:
            # self.debug("Counting permutations:")
            al = self.log.action("Counting permutations")
            ex_classes = self.exchangeable_classes(vars)
            disjoint_classes = self.disjoint(ex_classes)
            if len(ex_classes) == 1:
                # self.debug(f"Exchangeable variables...")
                self.log.detail(al, f"Exchangeable variables")
                sol = self.count_permutation_exchangeable(vars)
            # elif disjoint_classes:
            #     self.log("Different classes but disjoint...")
            #     sol = self.count_sequence_cartesian(vars)
            else:
                # self.debug("Splitting injectivity...")
                self.log.detail(al, "Splitting injectivity")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

    def count_sequence_any(self, var_list=None):
        """
        When variables are independent, solutions are from the cartesian product of the domain.
        Indistinguishable elements collapse to one in terms of possible choiches
        """
        vars = var_list if var_list is not None else self.vars
        count = One()
        al = self.log.action("Counting regular sequence")
        for v in vars:
            # the available choices are all dist. elems + each distinguishable property
            indist = v.elements.find(False)
            indist_sizes = self.get_sizes_indistinguishable(indist)
            indist_size = len([s for s in indist_sizes if s > 0])
            dist_size = v.size() - sum(
                indist_sizes
            )  # do not count each indistinguishable as a different choice
            # self.debug(f"\t{dist_size} different choices for {v}")
            self.log.detail(al, f"\t{dist_size} different choices for {v}")
            # but each set of indistinguishable counts as one possible choice
            choices = dist_size + indist_size
            count *= Count(choices, subproblems=0)
        # self.debug(f"\tDomain product: {count}")
        self.log.detail(al, f"\tDomain product: {count}")
        return count

    def count_permutation_exchangeable(self, var_list=None):
        """
        When all elements belong to the same set in a permutation, apply falling factorial
        However if sets have indistinguishable elements, account for it
        """
        vars = var_list if var_list is not None else self.vars
        n = len(vars)
        dom = vars[0]
        size = dom.size()
        extra = size - n
        al = self.log.action("Regular permutation")
        if extra >= 0:  # are there enough elements to fill the permutation?
            indist = dom.elements.find(False)
            dist = dom.elements.find(True)
            dist_lab = self.universe.get_label(dist)
            if indist.empty:
                self.log.detail(
                    al, f"Without indistinguishable: applying falling factorial"
                )
                if extra == 0:
                    count = Factorial(size, "Nr. orders for all objects")
                else:
                    num = Factorial(size, "Nr. orders for all objects")
                    denom = Factorial(extra, "Nr. excluded objects")
                    count = Divide(num, denom)
                count.set_histogram(self.histogram())
            else:
                self.log.detail(
                    al,
                    f"With indistinguishable elements: considering distributions of indistinguishable",
                )
                indist_sizes = self.get_sizes_indistinguishable(indist)
                dist_size = size - sum(indist_sizes)
                choices = self.get_group_choices([dist_size] + indist_sizes, n)
                labels = [self.universe.get_label(i) for i in indist]
                count = Zero()
                for choice in choices:
                    # each combination of choices is a valid disjoint subproblem
                    dist_choice = choice[0]
                    n_dist_choices = Binomial(
                        dist_size,
                        dist_choice,
                        f"Choose {dist_choice} of {dist_size} (distinguishable) {dist_lab} for {n} object(s)",
                    )
                    arrange_all = Factorial(
                        n,
                        f"Permutations of {n} {dom.name}",
                    )  # count = n! / n1!*n2!*...*nm!
                    arrange_indist_choices = [
                        Factorial(
                            ic,
                            f"Extra permutations of (indist.) {labels[i]}",
                        )
                        for i, ic in enumerate(choice[1:])
                    ]
                    prod_indist_choices = One()
                    for aic in arrange_indist_choices:
                        prod_indist_choices *= aic
                    count += Divide(
                        n_dist_choices * arrange_all,
                        prod_indist_choices,
                        histogram=self.histogram(),
                    )

        else:  # not enough elements in dom to make a n_vars long permutation
            count = Zero()
            self.log.detail(al, f"{len(vars)} different vars for {size} values!")
        return count

    def count_subsets(self, var_list=None):
        """
        Look at a subset problem and figure out what to do:
        - Subset with repetition: counting rule
        - Subset with alldiff (no repetition):
            - Exchangeable: counting rule (binomial coefficient)
            - Non-exchangeable: split/shatter
        """
        # self.log.action("Counting (multi)subsets")
        # self.log.detail(al,f"There are {self.fixed_choices} fixed elements")
        vars = var_list if var_list is not None else self.vars
        f = self.fixed_choices
        if self.config.type == "multisubset":
            ex_classes = self.exchangeable_classes(vars)
            if len(ex_classes) == 1:
                sol = self.count_multisubsets_exchangeable(vars[f:])
            else:
                sol = self.count_multisubsets_non_exchangeable(vars[f:])
        else:
            ex_classes = self.exchangeable_classes(vars)
            if len(ex_classes) == 1:
                sol = self.count_subsets_exchangeable(vars[f:])
            else:
                # self.debug("Splitting injectivity...")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

    def count_subsets_exchangeable(self, var_list=None):
        """
        When all elements belong to the same set in a permutation, apply binomial coefficient
        However if sets have indistinguishable elements, account for it
        """
        vars = var_list if var_list is not None else self.vars
        n = len(vars)
        dom = vars[0]
        size = dom.size()
        extra = size - n
        al = self.log.action("Counting regular subsets")
        if extra >= 0:
            indist = dom.elements.find(False)
            if indist.empty:
                count = Binomial(size, n, f"Choose {n} from {size} {dom.name} objects")
                self.log.detail(al, f"With binomial coefficient: ({size} {n})")
            else:
                indist_sizes = self.get_sizes_indistinguishable(indist)
                dist_size = size - sum(indist_sizes)
                choices = self.get_group_choices([dist_size] + indist_sizes, n)
                count = Zero()
                self.log.detail(
                    al, f"With indistinguishable elements: accounting for histograms"
                )
                self.log.detail(al, f"Indistinguishable groups: {indist}")
                for choice in choices:
                    self.log.detail(al, str(choice))
                    dist_choice = choice[0]
                    n_dist_choices = Binomial(
                        dist_size,
                        dist_choice,
                        f"Choose {dist_choice} from {dist_size} {dom.name} objects",
                    )
                    count += n_dist_choices
                # self.debug(f"Subsets with indistinguishable elements: {count}")
        else:  # not enough elements in dom to make a n sized subset
            self.log.detail(al, f"{len(vars)} different vars for {size} values!")
            count = Zero()
            # self.debug(f"{len(vars)} different vars for {size} values!")
        count.set_histogram(self.histogram())
        return count

    def get_sizes_indistinguishable(self, indistinguishable):
        """
        Given an interval containing the indistinguishable elements, returns
        the count of each group
        """
        indist_sizes = []
        for left, lower, upper, right in Int.to_data(indistinguishable):
            if not left:
                lower += 1
            if not right:
                upper -= 1
            if upper > lower:
                size = upper - lower + 1
            else:
                size = 0
            indist_sizes.append(size)
        return indist_sizes

    def get_group_choices(self, sizes, length):
        """
        For a sequence of given length, fix any way of picking n elements from each distinguishable property,
        a distinguishable property is either a distinguishable element (1 choice from it) or a set of
        indistinguishable elements (#set possible choices)
        sizes: size of each group of elements in the domain:
               first one is the group of distinguishable, others are groups of indistinguishable
        length: length of the sequence
        """
        if len(sizes) > 0:
            group_size = sizes[0]
            upper = min(length, group_size)
            lower = max(0, length - sum(sizes[1:]))
            indist_choices = []
            for s in range(
                lower, upper + 1
            ):  # for each set we can pick from 0 to #set elements
                rest_indist_choices = self.get_group_choices(sizes[1:], length - s)
                indist_choices += [[s] + ric for ric in rest_indist_choices]
            return indist_choices
        else:
            return [[]]

    def split_inj(
        self, split_class_vars, rest_classes_vars, split_class_ccs, rest_classes_ccs
    ):
        """
        Splits the problem under injectivity: we split between one split class and rest: split class has n variables and rest has some other classes/properties: we need to count how many "relevant" properties can end up in the split class. The relevant properties are the formulas that we have in rest: if we split between french and dutch, when splitting on french we consider the cases where n french speakers were dutch speakers, and adjust the domain of dutch in the "rest" problem.
        """
        split_class = split_class_vars[0]
        # self.debug(f"Split class: {split_class}")
        al = self.log.action(f"Splitting injectivity on class {split_class}")
        rest_domains = self.exchangeable_classes(rest_classes_vars).keys()
        relevant_rest_domains = [
            d for d in rest_domains if not split_class.disjoint(d)
        ]  # optimization: filter out disjoint domains
        indist_domains = [
            SetFormula(
                split_class.formula,
                split_class.elements[dom],
                split_class.universe,
                split_class.name,
            )
            for dom in split_class.elements.find(False)
        ]
        cases = self.relevant_cases(split_class, indist_domains + relevant_rest_domains)
        c = len(cases)
        if c == 0:
            # optimization: if all disjoint then split class subproblem already independent
            # self.debug(
            #     "Split class is disjoint from others: no need for injectivity split..."
            # )
            self.log.detail(
                al,
                "Split class is disjoint from others: no need for injectivity split",
            )
            count = self.split(
                split_class_vars, rest_classes_vars, split_class_ccs, rest_classes_ccs
            )
        else:
            n = len(split_class_vars)
            possible_partitions = self.integer_k_partitions(n, c)
            count = Zero()
            for partition in possible_partitions:
                possible_composition = set(itertools.permutations(partition))
                for id, composition in enumerate(possible_composition):
                    count += self.split_inj_case(
                        id,
                        composition,
                        cases,
                        split_class_vars,
                        rest_classes_vars,
                        split_class_ccs,
                        rest_classes_ccs,
                    )
        return count

    def split_inj_case(
        self,
        id,
        composition,
        cases,
        split_class_vars,
        rest_classes_vars,
        split_class_ccs,
        rest_classes_ccs,
    ):
        # self.debug(f"Case {n_case} of {cases}")
        split_count_formulas = split_class_ccs.copy()
        split_inj_formulas = [
            CCounting(cases[i], Int.singleton(composition[i]))
            for i in range(0, len(composition))
        ]
        split_count_formulas = split_count_formulas + split_inj_formulas
        # self.debug(
        #     f"Case {n} {split_class} are s.t. {split_count_formulas}:"
        # )
        try:
            # log info
            n = len(split_class_vars)
            v = split_class_vars[0]
            split_count_formulas = self.compact_ccs(split_count_formulas)
            description = f"Left split: case {n} {v} are s.t. {split_count_formulas}"

            split_class_count = self.solve_subproblem(
                split_class_vars,
                self.config,
                [],
                split_count_formulas,
                op="split-left",
                caption=description,
                id=f"{id+1}",
                shatter_id=id + 1,
            )
            if split_class_count.is_zero():
                # optimization: do not solve rest if split already unsat
                return Zero()
            for j, count_hist in enumerate(split_class_count.histograms):
                # consider all possible histograms for a left solution
                # i.e. compositions over relevant sets in a valid left configuration

                # update universe with left problem's histogram info
                c, hst = count_hist
                filtered = self.filter_domains(hst, rest_classes_vars)
                new_univ = self.filter_domains(hst, [self.universe])[0]

                # log info
                desc_hist = ", ".join([f"{n} {set}" for set, n in hst.items()])
                description = f"Right split removing {desc_hist}"

                # multiplication rule
                left_count = Count(c, histogram=hst)
                right_count = self.solve_subproblem(
                    filtered,
                    self.config,
                    [],
                    rest_classes_ccs,
                    universe=new_univ,
                    op="split-right",
                    caption=description,
                    id=f"{id+1}.H{j+1}",
                    shatter_id=id + 1,
                )
                # self.debug(
                #     f"Split result = {c} * {rest_count.count} = {partial_count * rest_count}"
                # )
                return left_count * right_count
        except Unsatisfiable as u:
            return Zero(tip=u.value)
            # self.debug(f"\t is unsat.")

    #########################
    ##  Level 2 methods ##
    #########################

    def count_exchangeable_class_partitions(
        self, size, n_elems, n_partitions, formula="universe"
    ):
        # https://math.stackexchange.com/questions/640558/how-many-ways-can-n-elements-be-partitioned-into-subsets-of-size-k
        self.log.action(
            f"Dividing {n_elems} from {formula} in {n_partitions} {self.config.type} of size {size}"
        )
        if size == 0:
            return One()
        pick = size * n_partitions
        count_pick = Factorial(
            pick,
            f"Nr. of ways to write a sequence of the {pick} picked objects where each block of {size} objects is a part",
        )
        # count_pick = Solution(
        #     math.factorial(pick),
        #     [],
        #     self.log,
        #     0,
        #     f"\\texttip{{ {pick}! }}{{Nr. of ways to write a sequence of the {pick} picked objects where each block of {size} objects is a part}}",
        # )
        count_pick_choice = Binomial(
            n_elems, pick, f"Pick {pick} of {formula} out of {n_elems}"
        )
        # pick_choice = math.comb(n_elems, pick)
        # count_pick_choice = Solution(
        #     pick_choice,
        #     [],
        #     self.log,
        #     0,
        #     f"\\binom{{ \\texttip{{ {n_elems} }}{{ Total nr. of objects}} }} {{ \\texttip{{ {pick} }}{{Nr. of objects to distribute }} }}",
        # )
        # counting rule
        # d1 = math.factorial(size) ** n_partitions
        # d2 = math.factorial(n_partitions)
        # count_d1 = Solution(
        #     d1,
        #     [],
        #     self.log,
        #     0,
        #     f"\texttip{{ ({size}!)^{n_partitions} }}{{ Each of the {n_partitions} has {size} internal orderings}}",
        # )
        # count_d2 = Solution(
        #     d2,
        #     [],
        #     self.log,
        #     0,
        #     f"\\texttip{{ {n_partitions}! }}{{ Nr. of external orderings of the {n_partitions} parts }}",
        # )
        d1 = Factorial(size, f"Internal ordering of each partition") ** Count(
            n_partitions, tip="Number of parts"
        )
        d2 = Factorial(
            n_partitions,
            tip=f"Nr. of external orderings of the {n_partitions} parts",
        )
        if self.config.type == "composition":
            count = Divide(
                count_pick, d1, tip="Total permutations divided by nr. internal orders"
            )
        else:
            count = Divide(
                count_pick,
                (d1 * d2),
                tip="Total permutations divided by nr. internal and external orders",
            )
        # TODO: with indist elements?
        # if self.type == "composition":
        #     count = count_pick // count_d1
        # else:
        #     count = count_pick // (count_d1 * count_d2)
        return count_pick_choice * count

    def count_fixed_partitions(
        self, relevant, var_list=None, caption="", id="1", log=None
    ):
        """
        If everything is fixed, we can count subsets like fixed sizes, i.e.
        for each partition pick the given number of p stated by (#p = n) and remove them from universe
        Account for exchangeability of variables in the process: pick at exchangeable class level
        """
        vars = var_list if var_list is not None else self.vars
        log = self.log if log is None else log
        al = log.action("Counting fixed partitions from histograms")
        for k in vars[0].histogram.keys():
            log.detail(al, str(k))
            log.add_relevant_set(k)
        choices = {rvs: rvs.size() for rvs in relevant}
        count = One()
        for i, v in enumerate(vars):
            al2 = log.action(
                f"Histogram for set {i+1} = {[val for val in v.histogram.values()]}"
            )
            for rvs in v.histogram:
                if not rvs.all_indistinguishable():
                    n = choices.get(rvs, 1)
                    k = v.histogram[rvs]
                    ec_choices = Binomial(n, k, tip=f"Choose {k} out of {n} {rvs}")
                    if rvs in choices:
                        choices[rvs] = n - k
                    count *= ec_choices
            log.detail(al2, f"Count = {count}")
        if self.config.type == "partition":
            part_histograms = [frozenset(v.histogram.items()) for v in vars]
            ex_hists = self.histogram(part_histograms)
            for hc in ex_hists:  # account for overcounting exchangeable choices (ugly)
                n_exchangeable = Factorial(ex_hists[hc], tip="Nr. exchangeable parts")
                count = Divide(count, n_exchangeable)
            log.detail(al, f"Count = {count}")
        return count

    def count_partitions(self, var_list=None, id=0):
        """
        Look at a partition problem and figure out what to do:
        - If there are no constraints and everything is exchangeable: counting rule (Stirling)
        - If there are size constraints solve with simpler solving method
        - If there are size constraints use the integer constraint solver method
        - If there are indistinguishable elements use the integer constraint solver method
        """
        vars = var_list if var_list is not None else self.vars
        p_log = ProblemLog(
            vars,
            universe=self.universe,
            caption=f"Constraint combination nr. {id}",
            id=id,
            level=self.lvl + 1,
            config=self.config,
            debug=self.debug,
        )
        n_elems = self.universe.size()
        ex_classes = self.exchangeable_classes(vars)
        counting_constraints = False
        for v in vars:
            counting_constraints = counting_constraints or len(v.ccs) > 0
        no_indist = self.universe.elements.find(False).empty
        if len(ex_classes) == 1:
            partition = vars[0]
            max_size = self.universe.size() - self.n_vars + 1
            unconstrained = CSize("any", Int.closed(1, max_size))
            if (
                not counting_constraints
                and partition.size == unconstrained
                and no_indist
            ):
                # self.debug(
                #     "Exchangeable variables and no constraints: Stirling numbers 2nd kind"
                # )
                al = p_log.action(
                    "Exchangeable variables and no constraints: Stirling numbers 2nd kind"
                )
                count = Stirling(
                    n_elems,
                    len(vars),
                    tip=f"Stirling number of the second kind: nr. partitions of {n_elems} distinguishable objecsts in {len(vars)} non-empty groups",
                )
                if self.config.type == "composition":
                    perm = Factorial(
                        len(vars), f"Nr. permutations of {len(vars)} groups"
                    )
                    count *= perm
            else:
                count = self.count_constrained_partitions(vars, log=p_log)
        else:
            count = self.count_constrained_partitions(vars, log=p_log)
        # self.debug(f"Count: {count}")
        return Solution(count, p_log)

    def count_constrained_partitions(self, var_list=None, log=None):
        """
        If all relevant quantities are fixed then count.
        Otherwise decompose intervals in constraints into equalities (accounting for all possible combinations)
        """
        # self.debug("There are constraints: consider relevant sets")
        log = self.log if log is None else log
        al = log.action("There are constraints: consider relevant sets")
        vars = var_list if var_list is not None else self.vars
        cases = self.universe.indistinguishable_subsets()
        for v in vars:
            cases = cases.union(v.relevant())
        cases = list(cases)
        if len(cases) == 1 and cases[0] == self.universe:
            relevant = cases
        else:
            relevant = self.relevant_cases_intersection(self.universe, cases)
        # self.debug(f"Relevant parts: {relevant}")
        log.detail(al, relevant)
        if self.fixed_vars(vars, relevant):
            return self.count_fixed_partitions(relevant, vars, log)
        else:
            return self.count_unfixed_partitions(relevant, vars, log)

    def fixed_vars(self, vars, relevant):
        """Check if everything relevant is fixed in each part/variable

        Args:
            vars ([LiftedSet]): parts
            relevant ([MultiSet]): list of relevant sets

        Return:
            Bool : all fixed
        """
        fixed = True
        for v in vars:
            fixed_size = is_singleton(v.size.values)
            fixed = fixed and fixed_size
            v.histogram = {rv_set: -1 for rv_set in relevant}
            if len(v.histogram) == 1 and fixed_size:  # register fix from size
                v.histogram[self.universe] = v.size.values.lower
            else:  # register fix from constraints
                for rv_set in relevant:
                    fixed_rv_set = False
                    for cc in v.ccs:
                        single = is_singleton(cc.values)
                        if rv_set == cc.formula and single:
                            v.histogram[rv_set] = cc.values.lower
                            fixed_rv_set = True
                    fixed = fixed and fixed_rv_set
        return fixed

    def fix_exchangeable_partititons(self, rv_set, n, prev_e=1, var_list=None):
        vars = var_list if var_list is not None else self.vars
        k = len(vars)
        # self.debug(f"Fixing {n} of {rv_set} on {k} exchangeable")
        # self.debug(f" {vars[0]}")
        al = self.log.action(f"Fixing {n} of {rv_set} on {k} exchangeable sets")
        self.log.detail(al, f"Exchangeable class: {vars[0]}")
        ips = self.feasible_relevant_partitions(rv_set, n, vars)
        valid = []
        for int_partition in ips:
            prop_vars = []
            for i, v in enumerate(vars):
                new_v = v.copy()
                new_v.histogram[rv_set] = int_partition[i]
                prop_vars.append(new_v)
            if len(prop_vars) == k:
                e = (
                    self.n_int_perms(int_partition)
                    if self.config.type == "composition"
                    else 1
                )
                # self.debug(f"{e} exchangeable {int_partition}")
                self.log.detail(al, f"{e} exchangeable fix: {int_partition}")
                valid.append((e * prev_e, prop_vars))
        return valid

    def fix_non_exchangeable_partititons(self, rv_set, prev_e=1, var_list=None):
        # self.debug(f"Fixing {rv_set} entities on non-exchangeable partitions:")
        self.log.action(f"Fixing {rv_set} entities on non-exchangeable partitions")
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        ics = self.feasible_relevant_compositions(ex_classes, rv_set)
        keys = list(ex_classes.keys())
        dists = []
        for distribution in ics:
            prop_vars = []
            # self.debug(
            #     f"Distribution of {rv_set} for {len(ex_classes)} classes: {distribution}"
            # )
            self.log.action(
                f"Distribution of {rv_set} for {len(ex_classes)} classes: {distribution}"
            )
            for j, n in enumerate(distribution):
                class_vars = ex_classes[keys[j]]
                prop_cvars = [v.copy() for v in class_vars]
                valid = self.fix_exchangeable_partititons(rv_set, n, 1, prop_cvars)
                if len(valid) > 0:
                    prop_vars.append(valid)
            if len(prop_vars) == len(ex_classes):  # otherwise something was unsat
                cases = list(self.product(prev_e, *prop_vars))
                dists += cases
        return dists

    def product(self, prev_e, *args):
        result = [(prev_e, [])]
        for pool in args:
            result = [(ex * ey, x + y) for ex, x in result for ey, y in pool]
        for prod in result:
            yield prod

    def fix_partitions(self, fixed, relevant):
        # fixed::(exchangeability constant, vars)
        if len(relevant) == 0:
            return [fixed]  # no more choices to account for
        else:
            al = self.log.action(f"Fixing {relevant} relevant sets")
            props = []
            e, vars = fixed
            rv_set = relevant[0]  # propagate one
            ex_classes = self.exchangeable_classes(vars)
            # self.debug(f"Propagating {rv_set} for {e}")
            self.log.detail(al, f"Propagating {rv_set} for set {e}")
            if len(ex_classes) == 1:
                prop_vars = self.fix_exchangeable_partititons(
                    rv_set, rv_set.size(), e, vars
                )
            else:
                prop_vars = self.fix_non_exchangeable_partititons(rv_set, e, vars)
            # self.debug(f"{rv_set} propagated")
            for pvs in prop_vars:
                prop_fix = self.fix_partitions(
                    pvs, relevant[1:]
                )  # propagate other relevants
                props += prop_fix
        return props

    def count_unfixed_partitions(self, relevant, var_list=None, log=None):
        vars = var_list if var_list is not None else self.vars
        log = self.log if log is None else log
        fixed = self.fix_partitions((1, vars), relevant)
        # self.debug(f"There are {len(fixed)} possible ways of fixing {relevant} domains")
        al = self.log.action("Counting unfixed partitions")
        self.log.detail(
            al, f"There are {len(fixed)} possible ways of fixing {relevant} domains"
        )
        count = Zero()
        for i, fix_vars in enumerate(fixed):
            e, prop_vars = fix_vars

            c_log = ProblemLog(
                prop_vars,
                universe=self.universe,
                caption=f"Fixing histograms, combination nr. {i+1}",
                id=i + 1,
                level=self.lvl + 1,
                config=self.config,
                debug=self.debug,
            )

            c_fix = self.count_fixed_partitions(
                relevant,
                prop_vars,
                caption=f"Fixed {i+1}",
                id=f"{self.log.id}.{str(i + 1)}",
                log=c_log,
            )
            c_log.count = c_fix.copy()
            log.add_subproblem("add", c_log)
            count += c_fix.with_choices(e)
        return count

    def feasible_relevant_partitions(self, rv_set, n, var_list=None):
        """
        Given a list of exchangerable variables LiftedSet, partitions n elements from rv_set.
        checking if there is some unsatisfiable constraint
        """
        vars = var_list if var_list is not None else self.vars
        int_partitions = self.integer_k_partitions(n, len(vars))
        fips = []
        al = self.log.action("Finding feasible relevant partitions")
        for int_partition in int_partitions:
            feasible = True
            for i, v in enumerate(vars):
                feasible = feasible and v.feasible(rv_set, int_partition[i])
            if feasible:
                fips.append(int_partition)
        # self.debug(f"{n} of {rv_set} have {len(fips)} feasible {len(vars)}-partitions")
        self.log.detail(
            al, f"{n} of {rv_set} have {len(fips)} feasible {len(vars)}-partitions"
        )
        return fips

    def feasible_relevant_compositions(self, ex_classes, rv_set):
        k = len(ex_classes)
        n = rv_set.size()
        int_partitions = self.integer_k_partitions(n, k)
        int_compositions = [
            list(c) for ip in int_partitions for c in set(itertools.permutations(ip))
        ]
        fics = []
        al = self.log.action("Finding feasible relevant compositions")
        for int_composition in int_compositions:
            feasible = True
            for i, v in enumerate(ex_classes.keys()):
                if feasible:
                    n_class = len(ex_classes[v])
                    f = v.feasible(rv_set, int_composition[i], n_class)
                    feasible = feasible and f
            if feasible:
                fics.append(int_composition)
        self.log.detail(
            al, f"{n} of {rv_set} have {len(fics)} feasible {k}-compositions"
        )
        return fics

    def n_int_perms(self, int_partition):
        # given the int partition returns the n!/n1!*n2!*...*nk! different permutations
        hst = self.histogram(int_partition)
        n = math.factorial(len(int_partition))
        p_classes = math.prod([math.factorial(hst[c]) for c in hst])
        return n // p_classes

    def propagate_ccs_partitions(self, vars, ccs, choices=1):
        """
        Given a list/set of variable partitions propagates counting constraints over partitions.
        If variables are exchangeable propagate, else shatter
        """
        # self.debug("Propagating size constraints...")
        al = self.log.action("Propagating size constraints")
        ex_classes = self.exchangeable_classes(vars)
        if len(ex_classes) == 1:  # pick a constraint and propagate
            cc = ccs[0]
            # self.debug(f"Variables are exchangeable: propagating {cc}")
            self.log.detail(al, f"Variables are exchangeable: propagating {cc}")
            if cc.values.upper == Int.inf:
                cc.values = cc.values.replace(upper=len(vars), right=Int.CLOSED)
            cases = []
            al2 = self.log.action("Shattering interval into single values")
            for i in Int.iterate(cc.values, step=1):
                self.log.detail(al2, f"Case  with {i} {cc.formula}")
                i_vars = vars.copy()
                cc_eq = CCounting(cc.formula, Int.singleton(i))
                try:
                    c_prop, p_vars = self.propagate_cc(
                        cc_eq, i_vars
                    )  # propagate one constraint
                    if len(ccs[1:]) > 0:
                        # self.debug("Shattering remaining...")
                        self.log.detail(al2, "Shattering remaining constraints")
                        shatter = self.propagate_ccs_partitions(
                            p_vars, ccs[1:], c_prop
                        )  # shatter others
                        cases += shatter
                    else:
                        cases += [(c_prop, p_vars)]
                except Unsatisfiable:
                    # self.debug("Propagation unsat")
                    self.log.detail(al2, "Propagation unsat")
            return cases
        else:  # shatter constraints
            # self.debug(f"Variables are not exchangeable: shattering size constraints")
            self.log.detail(
                al, f"Variables are not exchangeable: shattering size constraints"
            )
            scv, rcv = self.split_ex_classes(ex_classes)
            combs_split_class, combs_rest_classes = self.shatter_count_formulas(
                scv, rcv, ccs
            )
            vars_comb = []
            for i in range(0, len(combs_split_class)):
                cp_scv = scv.copy()
                cp_rcv = rcv.copy()
                prop_split_vars = self.propagate_ccs_partitions(
                    cp_scv, combs_split_class[i]
                )
                if len(prop_split_vars) == 1 and prop_split_vars[0][0] == 0:  # unsat
                    continue
                prop_rest_vars = self.propagate_ccs_partitions(
                    cp_rcv, combs_rest_classes[i]
                )
                if len(prop_rest_vars) == 1 and prop_rest_vars[0][0] == 0:  # unsat
                    continue
                ccs_combs = [
                    list(comb)
                    for comb in itertools.product(prop_split_vars, prop_rest_vars)
                ]
                for cc in ccs_combs:
                    s_problem, r_problem = cc
                    s_choices, svars = s_problem
                    r_choices, rvars = r_problem
                    c = choices * s_choices * r_choices
                    vs = svars + rvars
                    self.log.detail(al, f"Shattered as {vs}")
                    vars_comb.append((c, vs))
            return vars_comb

    def shatter_partitions(self):
        """
        Shatters partitions and combinations
        """
        if len(self.count_f) > 0:
            problems = self.propagate_ccs_partitions(self.vars, self.count_f)
        else:
            problems = [(1, vars)]
        count = Zero()
        # self.debug(
        #     f"Size constraints shattered, solving {len(problems)} combinations..."
        # )
        i = 1
        for c, p in problems:
            # self.debug(f"------ Combination {i} ------")
            i += 1
            prob_sol = self.count_partitions(p, i)
            self.log.add_subproblem("add", prob_sol.log)
            count += prob_sol.count.with_choices(c)
        return count

    def has_fixed_indist(self, var_list=None):
        vars = var_list if var_list is not None else self.vars
        ccs = vars[0].ccs
        indist_subsets = self.universe.indistinguishable_subsets()
        # indist = len(indist_subsets) > 0
        fixed = True
        for i in indist_subsets:
            fixed_set = False
            for cc in ccs:  # if a cc formula is a subset of i then i is partitioned
                fixed_set = fixed_set or (cc.formula in i)
            fixed = fixed and fixed_set  # if all sets are fixed then ok
        return fixed
