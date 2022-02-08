import portion as P
import math 
import itertools
from ortools.sat.python import cp_model
from portion.interval import singleton

from formulas import *
from configuration import Domain, LiftedSet, is_singleton

class Unsatisfiable(Exception):
    def __init__(self, value): 
        self.value = value 
  
    def __str__(self): 
        return(repr(self.value)) 


class Solution(object):
    """
    A solution is a set of pairs (count, histogram). 
    An histogram is used to keep track of used variables for injectivity constraint.
    For each histogram keep track of the corresponding number of (partial) solutions.
    """
    def __init__(self, count, histogram, subproblems =1):
        self.count = count
        self.histograms = [[count,histogram]]
        self.subproblems = subproblems
    
    def __add__(self, rhs):
        # if self.count == 0:
        #     return rhs
        # elif rhs.count == 0:
        #     return self
        # else:
        self.count += rhs.count
        self.histograms += rhs.histograms
        self.subproblems += rhs.subproblems
        return self

    def __mul__(self, rhs):
        histograms = []
        for c1, h1 in self.histograms:
            for c2, h2 in rhs.histograms:
                count = c1 * c2 
                histogram = self.add_histograms(h1, h2)
                histograms.append([count, histogram])
        prod = Solution(0,[])
        prod.count = self.count * rhs.count
        prod.histograms = histograms
        prod.subproblems = self.subproblems + rhs.subproblems
        return prod

    def __sub__(self, rhs):
        if rhs.count == 0:
            return self
        else:
            self.count -= rhs.count
            # self.histograms -= rhs.histograms
            self.subproblems += rhs.subproblems
            return self

    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return f"{self.count} ({self.subproblems} subproblems)"
    
    def add_histograms(self, h1, h2):
        sum = {}
        for f1 in h1:
            if f1 in h2:
                sum[f1] = h1[f1] + h2[f1]
            else:
                sum[f1] = h1[f1]
        for f2 in h2:
            if f2 not in h1:
                sum[f2] = h2[f2]
        return sum
    
    def with_choices(self, n_choices):
        updated_hst = []
        for c, hst in self.histograms:
            c = c*n_choices
            updated_hst.append([c,hst])
        self.count *= n_choices
        self.histograms = updated_hst
        return self


class SharpCSP(object):
    """
    Counts the possible assignments to a set of vars given choice/count constraints

    Attributes:
    vars : [DomainFormula]
        the variables representing the problem configuration
    type : str
        configuration type (same as Configuration)
    # agg_f : [AggFormula]
    #     aggrgate constraints on integers 
    choice_f : [PosFormula]
        choice formulas to enforce
    count_f : [CountFormulas]
        set of size constraints to enforce
    fixed_choices: int
        subsets only: number of variables already set by choice constraints
    lvl: int
        nesting level for logging 
    """

    def __init__(self, vars, type, choice_f, count_f, universe=None, lvl=0):
        # self.agg_f = agg_f
        self.choice_f = choice_f
        self.count_f = count_f
        self.fixed_choices = 0
        self.dolog = True
        self.lvl = lvl
        self.n_vars = len(vars)
        self.type = type
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
        if isinstance(chf, PosFormula):
            if chf.pos-1<self.n_vars:
                self.vars[chf.pos-1] = self.vars[chf.pos-1] & chf.dformula
                self.log("Choice set:", self.vars[chf.pos-1])
        # if isinstance(chf, InFormula):
        #     if self.fixed_choices > self.n_vars:
        #         raise Unsatisfiable("Too many elements for the subset size!")
        #     self.vars[self.fixed_choices] = chf.entity
        #     self.fixed_choices += 1
        #     self.log(f"Choice n.{self.fixed_choices} set:", self.vars[self.fixed_choices])

    def apply_count(self, cof, others):
        """
        Applies the constraints to satisfy 'current' count: for each admissible value n (lb<=n<=ub) 
        of the property p in cof, constrain the problem to n variables satisfying p and sum over the 
        different (independent) counts
        """
        self.log("Propagating ", cof)
        self.lvl += 1
        # count = Solution(0,[])
        # try:
        #     choices, prop_vars = self.propagate_cof(cof)
        #     self.vars = prop_vars
        #     count = self.split_on_constraints(choices, others)
        # except Unsatisfiable:
        #     pass
        # dom_size = cof.formula.size()+1 # n entities is the max count (open interval:+1)
        # if cof.values.upper == P.inf:
        #     # cof.values = cof.values.replace(upper = self.n_vars+1)
        #     cof.values = cof.values.replace(upper = dom_size)
        max_entities = min(self.n_vars,cof.formula.size())
        lb, ub = interval_closed(cof.values, ub_default=max_entities)
        all_values = P.closed(0,self.n_vars)
        out_values = all_values - cof.values
        count = Solution(0,[])
        if lb == ub:
            try:
                choices, prop_vars = self.propagate_cof(cof)
                self.vars = prop_vars
                count = self.split_on_constraints(choices, others)
            except Unsatisfiable:
                pass
        elif len(out_values)<len(cof.values):
            # ignore and then remove the unsat cases
            self.log(f"Relax {cof} and remove unsat")
            count_ignore = self.split_on_constraints(1, others)
            not_cof = CountingFormula(cof.formula, out_values)
            count_not = self.split_on_constraints(1, [not_cof] + others)
            count += count_ignore - count_not 
        else:
            # propositionalization of each i in the interval
            self.log(f"Expanding bounds {cof.values}...")
            values = cof.values.replace(upper=ub,right=P.CLOSED)
            for i in P.iterate(values, step = 1):
                cof_eq = CountingFormula(cof.formula, P.singleton(i))
                count_case = self.solve_subproblem(self.vars, self.type, {}, [cof_eq] + others)
                count += count_case
            self.log(f"Propagate lower bound #{cof.formula}={lb}...")
        self.lvl -=1
        return count
    
    def apply_counts(self):
        """
        If vars are exchangeable we can apply a count otherwise we need to split and consider the combinations of count constraints
        """
        ex_classes = self.exchangeable_classes()
        if self.type in  ["sequence", "permutation", "subset","multisubset"]:
            if len(ex_classes)>1:
                count = self.count_non_exchangeable(ex_classes)
            else:
                count = self.count_exchangeable()
        else:
            count = self.shatter_partitions()
        return count

    def choose_cof(self):
        """
        Heuristics to decide which counting constraint to propagate when variables are exchangeable
        """
        choice = self.count_f[0]
        others = self.count_f[1:]
        ex_class = self.vars[0]
        if isinstance(ex_class, DomainFormula):
            for cof in self.count_f:
                same_class = cof.formula == ex_class
                max = cof.values.lower == self.n_vars
                zero = cof.values.upper == 0 
                if same_class or max or zero:
                    choice = cof
                    others = self.count_f 
                    others.remove(cof)
                    return choice, others
            return choice, others
        else:
            return choice, others

    def compact_cofs(self, counts):
        """
        Simplify some counting constraints, i.e. turn (#p > n), (#p < m) into (#p in [n,m])
        """
        compact = []
        indexes = range(0,len(counts))
        remove = []
        updated = False
        for i in indexes :
            replaced = False
            cof1 = counts[i]
            if cof1.values.empty:
                raise Unsatisfiable("Empty set of valid counts")
            satisfied, not_satisfied, maybe = self.count_satisfied(cof1.formula)
            if len(maybe) == 0 and len(satisfied) in cof1.values: # remove cofs already satisfied
                    remove += [i]
            else:
                for j in [j for j in indexes if j>i]:
                    cof2 = counts[j]
                    if cof1 != cof2:
                        if cof1.formula == cof2.formula: # combination is not satisfiable
                            new_int = cof1.values & cof2.values
                            if new_int.empty:
                                raise Unsatisfiable("Conflicting constraints")
                            else:   # combination satisfiable: merge intervals into single cof
                                merged_cof = CountingFormula(cof1.formula, new_int)
                                compact.append(merged_cof)
                                remove += [i,j]
                    else: # cof1 == cof2:
                        #opposite formulas both 0
                        remove += [i]
        keep = [i for i in indexes if not i in remove]
        compact += list(map(lambda i: counts[i], keep))
        final =  len(remove) == 0
        result = compact if final else self.compact_cofs(compact)
        return result

    def count_exchangeable(self):
        """
        Take one count constraint and apply it
        """
        self.log("Counting exchangeable...")
        cof, others = self.choose_cof()
        count = self.apply_count(cof, others)
        return count
    
    def count_non_exchangeable(self, ex_classes):
        """
        If vars are not exchangeable we split the problem into class | rest
        Counting constraints need to be split according to the possible combinations 
        Example:
        Split classes french | students with cofs on dutch=2 and french=1
        for each cof prop==n the count can be split between the two classes with i elements on the left and n-i elements on the right, i.e. french==1 | french==0 and french==0 | french==1.
        But with many cofs we need to consider all the combinations of the two so:
        french==0 & dutch==0 | french==1 & dutch==2
        french==1 & dutch==0 | french==0 & dutch==2 ... 
        """
        split_class_vars, rest_classes_vars = self.split_ex_classes(ex_classes)
        combs_split_class, combs_rest_classes = self.shatter_count_formulas(split_class_vars, rest_classes_vars)
        tot_count = Solution(0,[])
        lvl = self.lvl
        for i in range(0,len(combs_rest_classes)):
            self.lvl = lvl
            try:
                comb_split_class = self.compact_cofs(combs_split_class[i])
                comb_rest_classes = self.compact_cofs(combs_rest_classes[i])
                self.log(f"Solving combination {i}: {comb_split_class} // {comb_rest_classes}")
                split_args = [split_class_vars, rest_classes_vars, list(comb_split_class), list(comb_rest_classes)]
                if self.type in ["permutation", "subset"]:
                    count = self.split_inj(*split_args)
                elif self.type in ["sequence", "multisubset"]:
                    count = self.split(*split_args)
                else:
                    count = self.split_partitions(*split_args)
                self.log(f"Split combination count: {count}")
            except Unsatisfiable:
                count = Solution(0,self.histogram())            
            tot_count += count
        self.lve = lvl
        self.log(f"Shatter count: {tot_count}")
        return tot_count

    def shatter_count_formulas(self, split_class_vars, rest_classes_vars, cofs=None):
        """
        Returns a shattering of the counting constraints for the given split.
        Given a split, for each size constraints compute the set of splits that shatter the constraint.
        Then combine each different cof shattering with the other constraints' possible splits (cartesian product).
        """
        self.log("Shattering size constraints...")
        count_f = self.count_f if cofs is None else cofs
        n_split = len(split_class_vars)
        n_rest = self.n_vars - n_split
        cofs_split_class = []
        cofs_rest_classes = []
        for cof in count_f:
            cases_split_class = []  # split for cof (left)
            cases_rest_classes = [] # split for cof (right)
            if cof.values.upper >= self.n_vars:
                interval_any = cof.values.replace(upper=n_split, right=P.CLOSED)    # adjust upper bound to number of vars
                any = CountingFormula(cof.formula, interval_any)    # add case >= n in split class then work on each =i <n
                cases_split_class.append(any)
                cof_rest_classes = CountingFormula(cof.formula, P.closed(0,n_rest))
                cases_rest_classes.append(cof_rest_classes)
                if cof.values.left == P.OPEN:
                    n_cases = cof.values.lower
                else:
                    n_cases = cof.values.lower-1
            else:
                n_cases = cof.values.upper
            for i in range(0,n_cases+1):
                cof_split_class = CountingFormula(cof.formula, P.singleton(i))
                cof_rest_classes = CountingFormula(cof.formula, cof.complement(i, n_rest))
                if self.is_feasible_split(split_class_vars[0],n_split,n_rest,cof_split_class,cof_rest_classes):
                    cases_split_class.append(cof_split_class)
                    cases_rest_classes.append(cof_rest_classes)
            cofs_split_class.append(cases_split_class)
            cofs_rest_classes.append(cases_rest_classes)
        combs_split_class =  list(itertools.product(*cofs_split_class))
        combs_rest_classes = list(itertools.product(*cofs_rest_classes))
        return combs_split_class, combs_rest_classes

    def count(self, vars=None):
        if self.type in ["sequence", "permutation"]:
            count = self.count_sequence(vars)
        elif self.type in ["subset", "multisubset"]:
            count = self.count_subsets(vars)
        elif self.type in ["partition","composition"]:
            count = self.count_partitions(vars)
        else:
            count = -1
        return count

    def get_lbub(self, interval):
        lower = interval.lower
        upper = interval.upper
        if interval.left == P.OPEN:
            lower +=1
        if interval.right == P.OPEN:
            upper -=1
        return lower, upper
    
    def count_satisfied(self, property, var_list = None):
        vars = var_list if var_list is not None else self.vars
        sat = []
        not_sat = []
        maybe = []
        if isinstance(property, DomainFormula):
            for i, v in enumerate(vars):
                if v in property:
                    sat.append(i)
                elif v.disjoint(property):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        elif isinstance(property, SizeFormula):
            for i, v in enumerate(vars):
                if v.size in property:
                    sat.append(i)
                elif v.size & property == SizeFormula("", P.empty()):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        else: 
            for i, v in enumerate(vars):
                satisfies = v.satisfies(property)
                if satisfies is None:
                    maybe.append(i)
                elif satisfies:
                    sat.append(i)
                else:
                    not_sat.append(i)
        return (sat, not_sat, maybe)

    def disjoint(self, classes):
        disj = True
        indexes = range(0,len(classes))
        for c1 in indexes:
            for c2 in [j for j in indexes if j>c1]:
                disj = disj and self.vars[c1].disjoint(self.vars[c2])
        return disj

    def exchangeable_classes(self, var_list = None):
        """
        Groups exchangeable variables in a dict
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
            if exclude.size()>0:
                df = dformula - exclude
        return df

    def filter_domains(self, cases, variables):
        self.log("Filtering domains...")
        dfs = variables
        for case in cases:
            self.log(f"  Case {case}")
            n = cases[case]
            if n> 0:
                pool = case
                for df in dfs[1:]:
                    if not df.disjoint(case):
                        pool = pool & df
                exclude = pool.take(n)
                if exclude.size()>0:
                    updated_dfs = []
                    for df in dfs:
                        self.log(f"  Filtering {n} {case} out of {df}:")
                        df = df-exclude
                        updated_dfs.append(df)
                    dfs = updated_dfs
                else:
                    dfs = [exclude] * n
        return dfs

    def get_vars(self, indexes):
        return [self.vars[v] for v in indexes]

    def histogram(self, var_list = None):
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
                yield a[:k + 2]
                x += 1
                y -= 1
            a[k] = x + y
            y = x + y - 1
            yield a[:k + 1]
    
    def integer_k_partitions(self, n, k):
        """
        Partitions n in k (possibly 0) integers
        """
        parts = []
        for p in self.integer_partitions(n):
            if len(p)<k:
                padded = p + [0]*(k-len(p))
                parts.append(padded)
            elif len(p) == k:
                parts.append(p)
        return parts       

    def is_feasible_split(self, split_class_var, n_split, n_rest, scof, rcof):
        """
        Checks if we ask to observe more properties than available variables
        Or things that split_class_var cannot satisfy
        """
        if scof.values.empty or rcof.values.empty:
            return False
        slb, sub =  interval_closed(scof.values, ub_default=n_split)
        rlb, rub =  interval_closed(rcof.values, ub_default=n_rest)
        if slb > n_split:
            return False
        if rlb > n_rest:
            return False
        disjoint = split_class_var.disjoint(scof.formula)
        included = scof.formula in split_class_var
        if disjoint and slb > 0:
            return False
        if included and sub < n_split:
            return False
        return True
    
    def log(self, s, *args):
        if self.dolog:
            strargs = " ".join(list(map(str,args)))
            flat = str(s) + strargs
            ind = "\t"*self.lvl
            lines = flat.split("\n")
            indented = list(map(lambda l: ind+l, lines))
            final = "\n".join(indented)
            print(final)
    
    def propagate(self, var, property):
        if isinstance(property, DomainFormula):
            return var & property
        elif isinstance(property, SizeFormula):
            return var & LiftedSet(self.universe, property)
        elif isinstance(property, CountingFormula):
            return var.add_cof(property)
        else:
            raise Exception(f"unexpected property type {property}: {type(property)}")
    
    def propagate_cof(self, cof, var_list = None):
        """
        Given a set of exchangeable variables, propagate one size constraint.
        - Check how many vars already satisfy, can satisfy, cannot satisfy the property
        - If there is no var that can satisfy, either is already sat or not
        - If there are vars that can satisfy, set m vars to satisfy the property where m is how many more entities there need to be, if set the others to not satisfy.
        return the number of exchangeable choices the propagation requires 
        """
        vars = var_list if var_list is not None else self.vars
        lb, ub = interval_closed(cof.values, ub_default=len(vars))
        satisfied, not_satisfied, maybe = self.count_satisfied(cof.formula, vars)
        diff = lb - len(satisfied)
        var_maybe = [vars[v] for v in maybe]
        if diff == 0 and len(maybe) == 0:
            self.log(cof, " already satisfied")
            return 1, vars
        elif len(maybe) == 0:
            self.log(cof, " is unsat here")
            raise Unsatisfiable("Unsat!")
        else: 
            v = var_maybe[0]
            self.log(f"{len(maybe)} exchangeable constrainable vars: {v}")
            if self.type == "subset" or self.type=="partition":
                n_choices = 1
            else:
                n_choices = math.comb(len(maybe), abs(diff))
            sat_f = self.propagate(v, cof.formula)
            not_sat_f = self.propagate(v, cof.formula.neg())
            for i in maybe:
                if diff<0:
                    vars[i]  = not_sat_f
                    diff += 1
                elif diff>0:
                    vars[i] = sat_f
                    diff -= 1
                else:
                    vars[i]  = not_sat_f
            if diff != 0:
                raise Unsatisfiable("Unsat!")
            else:
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
            if dom_base.size()>0:
                cases.append(dom_base) 
        return cases

    def relevant_cases(self, split_df, rest_classes):
        """
        Computes the cases for elements of type split_df w.rt. rest_classes. i.e. if we split on french and we have other classes like dutch and italian, 
        we need to consider for alldiff the cases where some of the french are both/neither dutch or italian or one of the two
        """
        self.log("Computing case combinations of relevant classes:")
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
        self.log(f"\t{res}")
        return res
       
    def solve(self, log=True):
        self.dolog = log
        self.log(self)
        for c in self.choice_f:
            self.apply_choice(c)
        try:
            self.count_f = self.compact_cofs(self.count_f)
        except Unsatisfiable:
            return Solution(0,[])
        if len(self.count_f) > 0:
            count = self.apply_counts()
        # elif len(self.agg_f) > 0:
        #     count = self.apply_aggs()
        else:
            count = self.count()
        self.log("=========")
        self.log("tot:" + str(count))
        return count

    def split(self, split_class_vars, rest_classes_vars, split_class_cofs, rest_classes_cofs):
        scv = [v.copy() for v in split_class_vars]
        rcv = [v.copy() for v in rest_classes_vars]
        self.log(f"Split class :")
        count = Solution(0, self.histogram())
        split_class_count = self.solve_subproblem(scv, self.type, [], split_class_cofs)
        if split_class_count != 0:
            self.log("Rest class: ")
            rest_classes_count = self.solve_subproblem(rcv, self.type, [], rest_classes_cofs)
            count = split_class_count * rest_classes_count
            self.log("==========")
        return count

    def solve_subproblem(self, vars, type, choice_constr, count_constr):
        self.log(f"\tSubproblem ({type}):")
        vars = [v.copy() for v in vars]
        ex_classes = self.exchangeable_classes(var_list=vars)
        if len(ex_classes) == 1:
            u = vars[0]
        else:
            u = self.universe
        subproblem = SharpCSP(vars, type, choice_constr, count_constr, u, self.lvl+1)
        try:
            count = subproblem.solve(self.dolog)
            return count
        except Unsatisfiable:
            self.log("\t==========\n\tUnsat: 0")
            return Solution(0, self.histogram())

    def split_on_constraints(self, n_choices, others):
        self.log("Splitting on other constraints...")
        if len(others) == 0:
            self.log("... no other constraints")
            count = self.count()
        else:
            count = self.solve_subproblem(self.vars, self.type, [], others)
        return count.with_choices(n_choices)
        
    def split_ex_classes(self, ex_classes):
        """
        Groups the exchangeable classes as left | right by taking one as left and uniting the others
        don't pick universe as split class
        """
        split_class = None
        rest_classes = []
        if self.type in ["partition", "composition"]:
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

    def stirling(self, n, k):
        computed = {}
        def stirling_aux(n, k):
            key = str(n) + "," + str(k)
            if key in computed.keys():
                return computed[key]
            if n == k == 0:
                return 1
            if (n > 0 and k == 0) or (n == 0 and k > 0):
                return 0
            if n == k:
                return 1
            if k > n:
                return 0
            result = k * stirling_aux(n - 1, k) + stirling_aux(n - 1, k - 1)
            computed[key] = result
            return result
        return stirling_aux(n,k)

    #######################
    ## Level 1 methods ##
    #######################
    
    def count_multisubsets(self, var_list = None):
        """
        Select with replacement from a set: if no indistinguishable elements use counting rule
        Otherwise for each indistinguishable set choose how many of its elements to put in the multisubset
        """
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        n = len(vars)
        dom = vars[0]
        indist = dom.elements.find(False)
        indist_sizes = self.get_sizes_indistinguishable(indist)
        indist_size = sum(indist_sizes)
        dist_size = dom.size() - indist_size + len(indist_sizes)
        # if indist_size == 0:
        if len(ex_classes) == 1:
            count = math.comb(dist_size+n-1,n)
        else:
            count = 0
            for exc in ex_classes:
                s = self.count_multisubsets(ex_classes[exc])
                count *= s.count
        sol = Solution(count, self.histogram())
        # else:
        #     choices = self.get_group_choices([dist_size] + indist_sizes, n) 
        #     count = 0
        #     for choice in choices:
        #         dist_choice = choice[0]
        #         n_dist_choices = math.comb(dist_size+dist_choice-1, dist_choice)
        #         # arrange_indist_choices = [math.factorial(ic) for ic in choice[1:]]
        #         count += n_dist_choices
        #     sol = Solution(count, self.histogram())
        return sol

    def count_sequence(self, var_list=None):
        """
        Look at a sequence problem and figure out what to do:
        - Sequence: cartesian product
        - Sequence with alldiff (permutation): 
            - Exchangeable: counting rule (factorial)
            - Non-exchangeable: split/shatter
        """
        vars = var_list if var_list is not None else self.vars
        if self.type == "sequence":
            self.log("Counting sequences:")
            sol = self.count_sequence_any(vars)
        else:
            self.log("Counting permutations:")
            ex_classes = self.exchangeable_classes(vars)
            disjoint_classes = self.disjoint(ex_classes)
            if len(ex_classes) == 1:
                self.log(f"Exchangeable variables...")
                sol = self.count_permutation_exchangeable(vars)
            # elif disjoint_classes:
            #     self.log("Different classes but disjoint...")
            #     sol = self.count_sequence_cartesian(vars)
            else:
                self.log("Splitting injectivity...")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

    def count_sequence_any(self, var_list = None):
        """
        When variables are independent, solutions are from the cartesian product of the domain.
        Indistinguishable elements collapse to one in terms of possible choiches 
        """
        vars = var_list if var_list is not None else self.vars
        count = 1
        for v in vars: # the available choices are all dist. elems + each distinguishable property
            indist = v.elements.find(False)
            indist_sizes = self.get_sizes_indistinguishable(indist) 
            indist_size = len([s for s in indist_sizes if s>0])
            dist_size = v.size() - sum(indist_sizes) # do not count each indistinguishable as a different choice
            self.log(f"\t{dist_size} different choices for {v}")
            count *= dist_size + indist_size # but each set of indistinguishable counts as one possible choice
        sol = Solution(count, self.histogram())
        self.log(f"\tDomain product: {count}")
        return sol

    def count_permutation_exchangeable(self, var_list = None):
        """
        When all elements belong to the same set in a permutation, apply falling factorial
        However if sets have indistinguishable elements, account for it
        """
        vars = var_list if var_list is not None else self.vars
        n = len(vars)
        dom = vars[0]
        size = dom.size()
        extra = size - n
        if extra >=0: # are there enough elements to fill the permutation?
            indist = dom.elements.find(False)
            if indist.empty:
                count = math.factorial(size) // math.factorial(extra)
                self.log(f"Falling factorial: {count}")
            else:
                indist_sizes = self.get_sizes_indistinguishable(indist)
                dist_size = size - sum(indist_sizes)
                choices = self.get_group_choices([dist_size] + indist_sizes, n)
                count = 0
                for choice in choices: # each combination of choices is a valid disjoint subproblem
                    dist_choice = choice[0]
                    n_dist_choices = math.comb(dist_size, dist_choice)
                    arrange_all = math.factorial(n) # count = n! / n1!*n2!*...*nm!
                    arrange_indist_choices = [math.factorial(ic) for ic in choice[1:]]
                    count += n_dist_choices * arrange_all // math.prod(arrange_indist_choices)
                self.log(f"Permutations with indistinguishable elements: {count}")
        else: # not enough elements in dom to make a n_vars long permutation
            count = 0
            self.log(f"{len(vars)} different vars for {size} values!")
        sol = Solution(count, self.histogram())
        return sol

    def count_subsets(self, var_list = None):
        """
        Look at a subset problem and figure out what to do:
        - Subset with repetition: counting rule 
        - Subset with alldiff (no repetition): 
            - Exchangeable: counting rule (binomial coefficient)
            - Non-exchangeable: split/shatter
        """
        self.log(f"There are {self.fixed_choices} fixed elements")
        vars = var_list if var_list is not None else self.vars
        f = self.fixed_choices
        if self.type == "multisubset":
            sol = self.count_multisubsets(vars[f:])
        else:
            ex_classes = self.exchangeable_classes(vars)
            if len(ex_classes) == 1:
                sol = self.count_subsets_exchangeable(vars[f:])
            else:
                self.log("Splitting injectivity...")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

    def count_subsets_exchangeable(self, var_list = None):
        """
        When all elements belong to the same set in a permutation, apply binomial coefficient
        However if sets have indistinguishable elements, account for it
        """
        vars = var_list if var_list is not None else self.vars
        n = len(vars)
        dom = vars[0]
        size = dom.size()
        extra = size - n
        if extra >=0:
            indist = dom.elements.find(False)
            if indist.empty:
                count = math.comb(size, n)
                self.log(f"Binomial coefficient: {count}")
            else:
                indist_sizes = self.get_sizes_indistinguishable(indist)
                dist_size = size - sum(indist_sizes)
                choices = self.get_group_choices([dist_size] + indist_sizes, n) 
                count = 0
                for choice in choices:
                    dist_choice = choice[0]
                    n_dist_choices = math.comb(dist_size, dist_choice)
                    count += n_dist_choices 
                self.log(f"Subsets with indistinguishable elements: {count}")
        else: # not enough elements in dom to make a n_vars sized subset
            count = 0
            self.log(f"{len(vars)} different vars for {size} values!")
        sol = Solution(count, self.histogram())
        return sol

    def get_sizes_indistinguishable(self, indistinguishable):
        """
        Given an interval containing the indistinguishable elements, returns 
        the count of each group
        """
        indist_sizes = []
        for left, lower, upper, right in P.to_data(indistinguishable):
            if not left:
                lower+= 1
            if not right:
                upper-= 1
            if upper > lower:
                size = upper-lower+1
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
            lower = max(0,length-sum(sizes[1:]))
            indist_choices = []
            for s in range(lower, upper+1): # for each set we can pick from 0 to #set elements
                rest_indist_choices = self.get_group_choices(sizes[1:], length-s)
                indist_choices += [[s] + ric for ric in rest_indist_choices] 
            return indist_choices
        else:
            return [[]]
                
    def split_inj(self, split_class_vars, rest_classes_vars, split_class_cofs, rest_classes_cofs):
        """
        Splits the problem under injectivity: we split between one split class and rest: split class has n variables and rest has some other classes/properties: we need to count how many "relevant" properties can end up in the split class. The relevant properties are the formulas that we have in rest: if we split between french and dutch, when splitting on french we consider the cases where n french speakers were dutch speakers, and adjust the domain of dutch in the "rest" problem.
        """
        self.lvl += 1
        split_class = split_class_vars[0]
        self.log(f"Split class: {split_class}")
        n = len(split_class_vars)
        rest_domains = list(self.exchangeable_classes(rest_classes_vars).keys())
        rest_domains = [d for d in rest_domains if not split_class.disjoint(d)] # optimization: filter out disjoint domains
        cases = self.relevant_cases(split_class, rest_domains)
        c = len(cases)
        if c == 0: # optimization: if all disjoint then split class subproblem already independent
            self.log("Split class is disjoint from others: no need for injectivity split...")
            count = self.split(split_class_vars, rest_classes_vars, split_class_cofs, rest_classes_cofs)
        else:
            count = Solution(0,[])
            for count_case in self.integer_k_partitions(n,c):
                ints = len(count_case)
                for n_case in set(itertools.permutations(count_case)):
                    self.log(f"Case {n_case} of {cases}")
                    split_count_formulas = split_class_cofs.copy()
                    split_inj_formulas = [CountingFormula(cases[i], P.singleton(n_case[i])) for i in range(0,len(n_case))]
                    split_count_formulas = split_count_formulas + split_inj_formulas
                    self.log(f"Case {n} {split_class} are s.t. {split_count_formulas}:")
                    try:
                        split_count_formulas = self.compact_cofs(split_count_formulas)
                        split_class_count = self.solve_subproblem(split_class_vars, self.type, [], split_count_formulas)
                        if split_class_count.count != 0: #optimization: do not solve rest if split already unsat
                            for c, hst in split_class_count.histograms:
                                partial_count = Solution(c, hst)
                                filtered = self.filter_domains(hst, rest_classes_vars)
                                rest_count = self.solve_subproblem(filtered, self.type, [], rest_classes_cofs)
                                self.log(f"Split result = {c} * {rest_count.count} = {partial_count * rest_count}")
                                count += partial_count * rest_count
                    except Unsatisfiable:
                        self.log(f"\t is unsat.")
        self.lvl -= 1
        return count

    #########################
    ##  Level 2 methods ##
    #########################

    def count_exchangeable_class_partitions(self, size, n_elems, n_partitions, formula = "universe"):
        # https://math.stackexchange.com/questions/640558/how-many-ways-can-n-elements-be-partitioned-into-subsets-of-size-k
        self.log(f"Dividing {n_elems} from {formula} in {n_partitions} {self.type} of size {size}")
        if size == 0:
            return 1
        pick = size*n_partitions
        pick_choice = math.comb(n_elems, pick)
        # counting rule
        d1 = math.factorial(size)**n_partitions
        d2 = math.factorial(n_partitions)
        count = math.factorial(pick)//(d1*d2)
        # se ho indist allora devo aggiustare 
        if self.type == "composition":
            count = count * math.factorial(n_partitions)
        return pick_choice*count
    
    def count_fixed_partitions(self, relevant, var_list = None):
        """
        If everything is fixed, we can count subsets like fixed sizes, i.e.
        for each partition pick the given number of p stated by (#p = n) and remove them from universe
        Account for exchangeability of variables in the process: pick at exchangeable class level
        """
        vars = var_list if var_list is not None else self.vars
        self.log("Counting fixed partitions on histogram ", [k for k in vars[0].histogram.keys()])
        choices = {rvs:rvs.size() for rvs in relevant}
        count = Solution(1,[])
        self.lvl += 1
        for i,v in enumerate(vars):
            self.log(f"Histogram V{i+1} = {[val for val in v.histogram.values()]}")
            for rvs in v.histogram:
                if not rvs.all_indistinguishable():
                    n = choices.get(rvs,1)
                    k = v.histogram[rvs]
                    ec_choices = math.comb(n,k) 
                    if rvs in choices: 
                        choices[rvs] = n - k  
                    count = count.with_choices(ec_choices)
        self.lvl -= 1
        if self.type=="partition":
            part_histograms = [frozenset(v.histogram.items()) for v in vars]
            ex_hists = self.histogram(part_histograms)
            for hc in ex_hists: # account for overcounting exchangeable choices (ugly)
                count.count = count.count // math.factorial(ex_hists[hc]) #should be int anyways
        self.log(f"\tcount = {count}")
        return count

    def count_partitions(self, var_list=None):
        """
        Look at a partition problem and figure out what to do:
        - If there are no constraints and everything is exchangeable: counting rule (Stirling)
        - If there are size constraints solve with simpler solving method
        - If there are size constraints use the integer constraint solver method
        - If there are indistinguishable elements use the integer constraint solver method
        """
        vars = var_list if var_list is not None else self.vars
        n_elems = self.universe.size()
        ex_classes = self.exchangeable_classes(vars)
        counting_constraints = False
        for v in vars:
            counting_constraints = counting_constraints or len(v.cofs)>0
        # se gli indistinguishable non sono fixed allora non posso andare su by_size
        if len(ex_classes) == 1: 
            partition = vars[0]
            max_size = self.universe.size() - self.n_vars + 1
            unconstrained = SizeFormula("any", P.closed(1, max_size))
            if not counting_constraints and partition.size == unconstrained:
                self.log("Exchangeable variables and no constraints: Stirling numbers 2nd kind")
                sol = self.stirling(n_elems, len(vars))
                if self.type == "composition":
                    sol *= math.factorial(len(vars))
                count = Solution(sol, self.histogram())
            else:
                count = self.count_constrained_partitions(vars)
        else:
            count = self.count_constrained_partitions(vars)
        self.log(f"Count: {count}")
        return count

    def count_constrained_partitions(self, var_list = None):
        """
        If all relevant quantities are fixed then count.
        Otherwise decompose intervals in constraints into equalities (accounting for all possible combinations)
        """
        self.log("There are constraints: consider relevant sets")
        vars = var_list if var_list is not None else self.vars
        cases = self.universe.indistinguishable_subsets()
        for v in vars:
            cases = cases.union(v.relevant())    
        cases = list(cases)
        if len(cases) == 1 and cases[0] == self.universe:
            relevant = cases
        else:
            relevant = self.relevant_cases_intersection(self.universe, cases)
        self.log(f"Relevant parts: {relevant}")
        fixed = True
        for v in vars: 
            fixed_size = is_singleton(v.size.values)
            fixed = fixed and fixed_size
            v.histogram = {rv_set: -1 for rv_set in relevant}
            if len(v.histogram) == 1 and fixed_size: # register fix from size
                v.histogram[self.universe] = v.size.values.lower
            else: # register fix from constraints
                for rv_set in relevant:
                    fixed_rv_set = False
                    for cof in v.cofs:
                        single = is_singleton(cof.values)
                        if rv_set == cof.formula and single:
                            v.histogram[rv_set] = cof.values.lower
                            fixed_rv_set = True
                    fixed = fixed and fixed_rv_set
        count = Solution(0,[])
        if fixed:
            count = self.count_fixed_partitions(relevant, vars)
        else:
            count = self.count_unfixed_partitions(relevant, vars)
        return count

    def fix_exchangeable_partititons(self, rv_set, n, prev_e=1,  var_list = None):
        vars = var_list if var_list is not None else self.vars
        k = len(vars)
        self.log(f"Fixing {n} of {rv_set} on {k} exchangeable")
        self.log(f" {vars[0]}")
        self.lvl += 1
        ips = self.feasible_relevant_partitions(rv_set, n, vars)
        valid = []
        for int_partition in ips:
            prop_vars = []
            for i, v in enumerate(vars):
                new_v = v.copy()
                new_v.histogram[rv_set] = int_partition[i]
                prop_vars.append(new_v)
            if len(prop_vars) == k:
                e = self.n_int_perms(int_partition) if self.type == "composition" else 1
                self.log(f"{e} exchangeable {int_partition}")
                valid.append((e*prev_e, prop_vars))
        self.lvl -= 1
        return valid

    def fix_non_exchangeable_partititons(self, rv_set, prev_e=1, var_list = None):
        self.log(f"Fixing {rv_set} entities on non-exchangeable partitions:")
        self.lvl += 1
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        ics = self.feasible_relevant_compositions(ex_classes, rv_set)
        keys = list(ex_classes.keys())
        dists = []
        for distribution in ics:
            prop_vars = []
            self.log(f"Distribution for {len(ex_classes)} classes: {distribution}")
            for j, n in enumerate(distribution):
                class_vars = ex_classes[keys[j]]
                prop_cvars = [v.copy() for v in class_vars]
                self.lvl += 1
                valid = self.fix_exchangeable_partititons(rv_set, n, prev_e, prop_cvars)  
                self.lvl -= 1
                if len(valid) > 0:                  
                    prop_vars.append(valid)
            if len(prop_vars) == len(ex_classes): # otherwise something was unsat
                cases = list(self.product(*prop_vars))
                dists += cases
        self.lvl -= 1
        return dists
    
    def product(self, *args):
        result = [(1,[])]
        for pool in args:
            result = [(ex*ey, x+y) for ex, x in result for ey, y in pool]
        for prod in result:
            yield prod

    def fix_partitions(self, fixed, relevant):
        # fixed::(exchangeability constant, vars)
        if len(relevant) == 0:
            return [fixed] # no more choices to account for
        else:
            props = []
            e, vars = fixed
            rv_set = relevant[0] # propagate one
            ex_classes = self.exchangeable_classes(vars)
            if len(ex_classes) == 1:
                prop_vars = self.fix_exchangeable_partititons(rv_set, rv_set.size(), e, vars)
            else: 
                prop_vars = self.fix_non_exchangeable_partititons(rv_set, e, vars)
            for pvs in prop_vars:
                prop_fix = self.fix_partitions(pvs, relevant[1:]) # propagate other relevants
                props += prop_fix
        return props
                
    def count_unfixed_partitions(self, relevant, var_list = None):
        vars = var_list if var_list is not None else self.vars
        fixed = self.fix_partitions( (1, vars), relevant)
        self.log(f"There are {len(fixed)} possible ways of fixing {relevant} domains")
        count = Solution(0,[])
        for fix_vars in fixed:
            e, prop_vars = fix_vars
            c_fix = self.count_fixed_partitions(relevant, prop_vars)
            count += c_fix.with_choices(e)
        return count

    def feasible_relevant_partitions(self, rv_set, n, var_list = None):
        """
        Given a list of exchangerable variables LiftedSet, partitions n elements from rv_set.
        checking if there is some unsatisfiable constraint 
        """
        vars = var_list if var_list is not None else self.vars
        int_partitions = self.integer_k_partitions(n, len(vars))
        fips = [] 
        for int_partition in int_partitions:
            feasible = True
            for i, v in enumerate(vars):
                feasible = feasible and v.feasible(rv_set, int_partition[i])
            if feasible:
                fips.append(int_partition)
        self.log(f"{n} of {rv_set} have {len(fips)} feasible {len(vars)}-partitions")
        return fips

    def feasible_relevant_compositions(self, ex_classes, rv_set):
        k = len(ex_classes)
        n = rv_set.size()
        int_partitions = self.integer_k_partitions(n, k)
        int_compositions = [list(c) for ip in int_partitions for c in set(itertools.permutations(ip))]
        fics = []
        for int_composition in int_compositions:
            feasible = True
            for i, v in enumerate(ex_classes.keys()):
                n_class = len(ex_classes[v])
                f = v.feasible(rv_set, int_composition[i], n_class)
                feasible = feasible and f
            if feasible:
                fics.append(int_composition)
        return fics

    def n_int_perms(self, int_partition):
        # given the int partition returns the n!/n1!*n2!*...*nk! different permutations
        hst = self.histogram(int_partition)
        n = math.factorial(len(int_partition))
        p_classes = math.prod([math.factorial(hst[c]) for c in hst])
        return n // p_classes

    def propagate_cofs_partitions(self, vars, cofs, choices=1):
        """
        Given a list/set of variable partitions propagates counting constraints over partitions.
        If variables are exchangeable propagate, else shatter
        """
        self.log("Propagating size constraints...")
        self.lvl += 1
        ex_classes = self.exchangeable_classes(vars)
        if len(ex_classes) == 1: # pick a constraint and propagate
            cof = cofs[0]
            self.log(f"Variables are exchangeable: propagating {cof}")
            if cof.values.upper == P.inf:
                cof.values=cof.values.replace(upper=len(vars), right = P.CLOSED)
            cases = []
            for i in P.iterate(cof.values, step = 1):
                i_vars = vars.copy()
                cof_eq = CountingFormula(cof.formula, P.singleton(i))
                # dl = self.dolog
                try:
                    # self.dolog = False
                    c_prop, p_vars = self.propagate_cof(cof_eq, i_vars) # propagate one constraint
                    if len(cofs[1:])>0:
                        self.log("Shattering remaining...")
                        shatter = self.propagate_cofs_partitions(p_vars, cofs[1:], c_prop) # shatter others
                        cases += shatter
                    else:
                        cases += [(c_prop, p_vars)]
                    # self.dolog = dl
                    # self.lvl -= 1
                except Unsatisfiable:
                    self.log("Propagation unsat")
                    # self.dolog = dl
                    self.lvl -= 1
            self.lvl -= 1
            return cases
        else: # shatter constraints
            self.log(f"Variables are not exchangeable: shattering size constraints")
            scv, rcv = self.split_ex_classes(ex_classes)
            combs_split_class, combs_rest_classes = self.shatter_count_formulas(scv, rcv, cofs)
            vars_comb = []
            for i in range(0,len(combs_split_class)):
                cp_scv = scv.copy()
                cp_rcv = rcv.copy()
                prop_split_vars = self.propagate_cofs_partitions(cp_scv, combs_split_class[i])
                if len(prop_split_vars) == 1 and prop_split_vars[0][0] == 0: #unsat
                    continue
                prop_rest_vars = self.propagate_cofs_partitions(cp_rcv, combs_rest_classes[i])
                if len(prop_rest_vars) == 1 and prop_rest_vars[0][0] == 0: #unsat
                    continue
                cofs_combs = [list(comb) for comb in itertools.product(prop_split_vars, prop_rest_vars)]
                for cc in cofs_combs:
                    s_problem, r_problem = cc
                    s_choices, svars = s_problem
                    r_choices, rvars = r_problem
                    c = choices * s_choices * r_choices
                    vars_comb.append((c, svars+rvars))
            self.lvl -= 1
            return vars_comb     

    def shatter_partitions(self):
        """
        Shatters partitions and combinations 
        """
        if len(self.count_f) > 0:
            problems = self.propagate_cofs_partitions(self.vars, self.count_f)
        else:
            problems = [(1, vars)]  
        count = Solution(0, [])
        self.log(f"Size constraints shattered, solving {len(problems)} combinations...")
        i = 1
        for c, p in problems:
            self.log(f"------ Combination {i} ------")
            i += 1
            self.lvl += 1
            count += self.count_partitions(p).with_choices(c)
            self.lvl -= 1
        return count

    def has_fixed_indist(self, var_list=None):
        vars = var_list if var_list is not None else self.vars
        cofs = vars[0].cofs
        indist_subsets = self.universe.indistinguishable_subsets()
        # indist = len(indist_subsets) > 0
        fixed = True
        for i in indist_subsets:
            fixed_set = False
            for cof in cofs: # if a cof formula is a subset of i then i is partitioned 
                fixed_set = fixed_set or (cof.formula in i)
            fixed = fixed and fixed_set # if all sets are fixed then ok
        return fixed
    

class SolutionCounter(cp_model.CpSolverSolutionCallback):
    """
    Count intermediate solutions.
    Interface to ortools solver to generate a valid combination of counting constraints with fixed sizes
    (all of the form #p=n)
    """

    def __init__(self, type, variables, cases, universe, lvl):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__type = type
        self.__variables = variables
        self.__cases = cases
        self.__universe = universe
        self.__lvl = lvl
        self.__solution_count = Solution(0, [])

    def OnSolutionCallback(self):
        var_partitions = []
        for i, partition in enumerate(self.__variables):
            s = sum([self.Value(v) for v in partition])
            size = SizeFormula(f"s{i}", P.singleton(s))
            ls = LiftedSet(f"p{i}", size)
            for case, var in enumerate(partition):
                v = self.Value(var)
                value = P.singleton(v)
                cof = CountingFormula(self.__cases[case], value)
                ls.cofs.append(cof)
            var_partitions.append(ls)
        # recursively solve the problem, now with fixed value constraints
        try:
            csp = SharpCSP(var_partitions, self.__type, [], [], self.__universe, self.__lvl +1) 
            sol = csp.solve(True)
            self.__solution_count += sol
        except Exception as e:
            print(e)

    def SolutionCount(self):
        return self.__solution_count

