import portion as P
import math 
import itertools
import operator
from ortools.sat.python import cp_model

from problog.logic import Constant

from formulas import *
from structure import Domain, LiftedSet


class Unsatisfiable(Exception):
    def __init__(self, value): 
        self.value = value 
  
    def __str__(self): 
        return(repr(self.value)) 

class SolutionCounter(cp_model.CpSolverSolutionCallback):
  """Count intermediate solutions."""

  def __init__(self, variables, cases):
    cp_model.CpSolverSolutionCallback.__init__(self)
    self.__variables = variables
    self.__cases = cases
    self.__solution_count = Solution(0, [])

  def OnSolutionCallback(self):
    sol = Solution(1,[])
    case_elems = [case.size() for case in self.__cases]
    for i, partition in enumerate(self.__variables):
        print(f"p{i}:", end = ' ')
        for case,v in enumerate(partition):
            val = self.Value(v)
            n = case_elems[case]
            print(f"{v} = {val}", end = ' ')
            choices = math.comb(n,val)
            case_elems[case] -= val
            sol = sol.with_choices(choices)
        print()
    print("\t",sol)
    self.__solution_count += sol

  def SolutionCount(self):
    return self.__solution_count

class Solver(object):
    """
    Sets up the decision variables from the problem
    """
    def __init__(self,problem):
        self.problem = problem
        self.universe = problem.domains[problem.universe]
        self.size = problem.structure.size
        self.type = problem.structure.type

    def solve(self, log=True):
        count = Solution(0,[])
        var_dom = DomainFormula(self.universe, Constant(self.problem.universe), self.universe)
        for n in self.size:
            if self.type in ["sequence", "subset"]:
                vars = [var_dom]*n
            else:
                ub_size = self.universe.size() - n + 1
                size = SizeFormula("universe", P.closed(1,ub_size))
                vars = [LiftedSet("liftedset", size, var_dom)]*n
            csp = SharpCSP(vars, self.type, self.problem.choice_formulas, self.problem.count_formulas, self.problem.structure.spec)
            count += csp.solve(log)
        return count

class Solution(object):
    """
    A solution is a set of pairs (count, histogram). 
    An histogram is used to keep track of used variables for injectivity constraint.
    For each histogram keep track of the corresponding number of (partial) solutions.
    """
    def __init__(self, count, histogram):
        self.count = count
        self.histograms = [[count,histogram]]
    
    def __add__(self, rhs):
        if self.count == 0:
            return rhs
        elif rhs.count == 0:
            return self
        else:
            self.count += rhs.count
            self.histograms += rhs.histograms
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
        return prod
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return str(self.count)
    
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
        the variables representing the problem structure
    type : str
        structure type (same as Structure)
    choice_f : [PosFormula/InFormula]
        choice formulas to enforce
    count_f : [CountFormulas]
        set of counting formulas to enforce
    alt_type : boolean
        same as Structure
    fixed_choices: int
        subsets only: number of variables already set by choice constraints
    lvl: int
        nesting level for logging 
    """

    def __init__(self, vars, type, choice_f, count_f, alt_type, lvl=0):
        self.alt_type = alt_type
        self.choice_f = choice_f
        self.count_f = count_f
        self.fixed_choices = 0
        self.dolog = True
        self.lvl = lvl
        self.n_vars = len(vars)
        self.type = type
        self.vars = vars

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
        str += "----------"
        return str

    def apply_choice(self, chf):
        """
        Sets variables according to the choice constraints
        """
        if isinstance(chf, PosFormula):
            self.vars[chf.pos-1] = self.vars[chf.pos-1] & chf.dformula
            self.log("Choice set:", self.vars[chf.pos-1])
        if isinstance(chf, InFormula):
            if self.fixed_choices > self.n_vars:
                raise Unsatisfiable("Too many elements for the subset size!")
            self.vars[self.fixed_choices] = chf.entity
            self.fixed_choices += 1
            self.log(f"Choice n.{self.fixed_choices} set:", self.vars[self.fixed_choices])

    def apply_count(self, cof, others):
        """
        Applies the constraints to satisfy 'current' count: for each admissible value n (lb<=n<=ub) 
        of the property p in cof constrain the problem to n variables satisfying p and sum over the 
        different (independent) counts
        """
        self.log("Propagating ", cof)
        self.lvl += 1
        lb = cof.values.lower
        ub = cof.values.upper
        if ub == P.inf:
            cof.values = cof.values.replace(upper = self.n_vars+1)
        count = Solution(0,[])
        if lb == ub:
            try:
                count = self.count_eq(cof, others)
            except Unsatisfiable:
                pass
        else:
            self.log(f"Expanding bounds {cof.values}...")
            for i in P.iterate(cof.values, step = 1):
                cof_eq = CountingFormula(cof.formula, P.singleton(i))
                count_case = self.solve_subproblem(self.vars, self.type, {}, [cof_eq] + others, self.alt_type)
                count += count_case
        self.lvl -=1
        return count
    
    def apply_counts(self):
        """
        If vars are exchangeable we can apply a count otherwise we need to split and consider the combinations of count constraints
        """
        ex_classes = self.exchangeable_classes()
        if len(ex_classes)>1:
            count = self.count_non_exchangeable(ex_classes)
        else:
            count = self.count_exchangeable()
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
        compact = []
        indexes = range(0,len(counts))
        remove = []
        updated = False
        for i in indexes :
            replaced = False
            cof1 = counts[i]
            satisfied, not_satisfied, maybe = self.count_satisfied(cof1.formula)
            if len(maybe) == 0 and len(satisfied) in cof1.values:
                    remove += [i]
            else:
                for j in [j for j in indexes if j>i]:
                    cof2 = counts[j]
                    if cof1 != cof2:
                        if cof1.formula == cof2.formula:
                            new_int = cof1.values & cof2.values
                            if new_int.empty:
                                raise Unsatisfiable("Conflicting constraints")
                            else:
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

    def count_eq(self, cof, others):
        """
        - Check how many vars already satisfy, can satisfy, cannot satisfy the property
        - If there is no var that can satisfy, either is already sat or not
        - If there are vars that can satisfy, set m vars to satisfy the property where m is how many more entities there need to be, set the others to not satisfy.
        - If we had enough variables go on with other constraints
        """
        satisfied, not_satisfied, maybe = self.count_satisfied(cof.formula)
        goal = cof.values.lower
        diff = goal - len(satisfied)
        var_maybe = self.get_vars(maybe)
        ex_classes = self.exchangeable_classes(var_maybe)
        if len(maybe) == 0:
            n_choices = 1
            if diff == 0:
                self.log(cof, " already satisfied")
                return self.split_on_constraints(n_choices, others)
            else:
                self.log(cof, " is unsat here")
                raise Unsatisfiable("Unsat!")
        else: #len(ex_classes) == 1:
            ex_class = next(iter(ex_classes))
            self.log(f"{len(maybe)} exchangeable constrainable vars: {ex_class}")
            if self.type == "subset" and self.alt_type == False or self.type=="partition":
                n_choices = 1
            else:
                n_choices = math.comb(len(maybe), abs(diff))

            sat_f = self.propagate(ex_class, cof.formula)
            not_sat_f = self.propagate(ex_class, cof.formula.neg())
            for i in maybe:
                if diff<0:
                    self.vars[i]  = not_sat_f
                    diff += 1
                elif diff>0:
                    self.vars[i] = sat_f
                    diff -= 1
                else:
                    self.vars[i]  = not_sat_f
            if diff != 0:
                raise Unsatisfiable("Unsat!")
            else:
                return self.split_on_constraints(n_choices, others)

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
        self.log("Counting non-exchangeable...")
        split_class_vars, rest_classes_vars = self.split_ex_classes(ex_classes)
        n_split = len(split_class_vars)
        n_rest = self.n_vars - n_split
        cofs_split_class = []
        cofs_rest_classes = []
        for cof in self.count_f:
            cases_split_class = []
            cases_rest_classes = []
            if cof.values.upper >= self.n_vars:
                interval_any = cof.values.replace(upper=n_split, right=P.CLOSED) # adjust upper bound to number of vars
                any = CountingFormula(cof.formula, interval_any) # add case >= n in split class then work on each =i <n
                cases_split_class.append(any)
                cof_rest_classes = CountingFormula(cof.formula, P.closed(0,n_rest))
                cases_rest_classes.append(cof_rest_classes)
                n_cases = cof.values.lower-1
            else:
                n_cases = cof.values.upper
            for i in range(0,n_cases+1):
                cof_split_class = CountingFormula(cof.formula, P.singleton(i))
                cof_rest_classes = CountingFormula(cof.formula, cof.complement(i, n_rest))
                if self.is_feasible_split(n_split,n_rest,cof_split_class,cof_rest_classes):
                    cases_split_class.append(cof_split_class)
                    cases_rest_classes.append(cof_rest_classes)
            cofs_split_class.append(cases_split_class)
            cofs_rest_classes.append(cases_rest_classes)
        combs_split_class =  list(itertools.product(*cofs_split_class))
        combs_rest_classes = list(itertools.product(*cofs_rest_classes))
        tot_count = Solution(0,[])
        lvl = self.lvl
        for i in range(0,len(combs_rest_classes)):
            self.lvl = lvl
            # comb_split_class = combs_split_class[i]
            # comb_rest_classes = combs_rest_classes[i]
            try:
                comb_split_class = self.compact_cofs(combs_split_class[i])
                comb_rest_classes = self.compact_cofs(combs_rest_classes[i])
                self.log(f"Solving combination {i}: {comb_split_class} // {comb_rest_classes}")
                split_args = [split_class_vars, rest_classes_vars, list(comb_split_class), list(comb_rest_classes)]
                if self.type in ["sequence", "subset"]:
                    if self.alt_type: # Fix true/false difference in language
                        count = self.split_inj(*split_args)
                    else:
                        count = self.split(*split_args)
                else:
                    count = self.split_partitions(*split_args)
                self.log(f"Split combination count: {count}")
            except Unsatisfiable:
                count = 0            
            tot_count += count
        self.lve = lvl
        self.log(f"Shatter count: {tot_count}")
        return tot_count

    def count(self, vars=None):
        if self.type=="sequence":
            count = self.count_sequence(vars)
        elif self.type=="subset":
            count = self.count_subsets(vars)
        elif self.type=="partition":
            count = self.count_partitions(vars)
        else:
            pass
        return count

    def count_partitions(self, var_list=None):
        self.log("Counting partitions:")
        vars = var_list if var_list is not None else self.vars
        ex_classes = self.exchangeable_classes(vars)
        counting_constraints = False
        for v in vars:
            counting_constraints = counting_constraints or len(v.constraints)>0
        if len(ex_classes) == 1: 
            partition = vars[0]
            max_size = partition.source.domain.size()-self.n_vars+1
            unconstrained = SizeFormula(1,max_size)
            if not counting_constraints and partition.size == unconstrained: 
                count = self.stirling(n, self.n_vars)
                sol = Solution(count, self.histogram())
            elif not counting_constraints and partition.size != unconstrained:
                count = self.count_partitions_by_size()
            else:
                count = self.count_ne_partitions(ex_classes, vars[0].source)
        elif counting_constraints:
            count = self.split_partitions()
        else: # sizes different but no extra constraint:
            count = self.count_partitions_by_size()
        return count

    def count_partitions_by_size(self, var_list = None, n = None):
        """
        Counts the number of partitions under the assumption that partitions do not have counting constraints.
        Fixes a size s of the subset, then recursively solves the problem on the other subsets where there are
        s elements less available to fill the partition. n keeps track of the available elements
        """
        vars = var_list if var_list is not None else self.vars
        v = vars[0]
        if n is None:
            n = v.source.size()
        if len(vars) > 1 and n>len(vars):
            var_count = Solution(0,[])
            for s in v.size:
                choices = math.comb(n,s)
                if choices > 0:
                    rest = self.count_partitions_by_size(vars[1:], n-s)
                var_count += rest.with_choices(choices)
            self.log(f"{var_count} solutions for var {v}")
        elif n==len(vars) and 1 in v.size or len(vars) == 1 and n in v.size:
            var_count = Solution(1, self.histogram())
        else:
            var_count = Solution(0,self.histogram())
        return var_count

    def count_ne_partitions(self, ex_classes, univ):
        n = univ.size()
        size_partitions = self.integer_k_partitions(n, self.n_vars)
        valid_size_partitions = [sp for sp in size_partitions if self.is_feasible_partition(sp)]
        relevant = set()
        for c in ex_classes:
            relevant = relevant.union(c.relevant())
        relevant = list(relevant)
        cases = self.relevant_cases(univ, relevant)
        n_cases = len(cases)
        count = Solution(0,[])
        for vsp in valid_size_partitions:
            p_model =  cp_model.CpModel()
            # decision variables
            case_vars = []
            for i, v in enumerate(self.vars):
                cvs = []
                for j,c in enumerate(cases):
                    cvs += [p_model.NewIntVar(0, c.size() ,f"c{i}{j}")]
                case_vars.append(cvs)
            for i,c in enumerate(cases):
                ith_case_vars = [cv[i] for cv in case_vars]
                p_model.Add(sum(ith_case_vars) == c.size())
            # local constraints
            for i,v in enumerate(self.vars):
                shatter_vars = case_vars[i]
                p_model.Add(sum(shatter_vars) == vsp[i])
                # counting constraints
                for cof in v.constraints:
                    subcases = [j for j,c in enumerate(cases) if c in cof.formula]
                    cof_vars = [shatter_vars[s] for s in subcases]
                    if cof.values.atomic:
                        lb, ub = self.get_lbub(cof.values)
                        p_model.Add(sum(cof_vars) >= lb)
                        p_model.Add(sum(cof_vars) <= ub)
                    else:
                        for vals_interval in cof.values:
                            lb, ub = self.get_lbub(vals_interval)
                            p_model.Add(sum(cof_vars) >= lb)
                            p_model.Add(sum(cof_vars) <= ub)

            solver = cp_model.CpSolver()
            solution_counter = SolutionCounter(case_vars, cases)
            status = solver.SearchForAllSolutions(p_model, solution_counter)
            count += solution_counter.SolutionCount()
        return count

    def get_lbub(self, interval):
        lower = interval.lower
        upper = interval.upper
        if interval.left == P.OPEN:
            lower +=1
        if interval.right == P.OPEN:
            upper +=1
        return lower, upper
    
    # def count_partitions_by_cofs(self, var_list = None):
    #     vars = var_list if var_list is not None else self.vars
    #     v = vars[0]
    #     relevant = list(v.relevant()) # exchangeable, altrimenti devo guardare anche le altre
    #     cofs = v.constraints
    #     splits = self.get_partition_splits(self, relvant)
    #     size_partitions = self.integer_k_partitions(n, len(vars))
    #     valid_size_partitions = [sp for sp in size_partitions if self.is_feasible_partition(v.size, sp)]
    #     for i, v in enumerate(self.vars):
    #         # stabilire gli histogram validi
    #         # ricorsione sul problema fissato l'histogram
    #         # filtrare dom dall'histogram fissato e continuare
    #     return Solution(0,[])

    # def get_partition_histograms(self, cases, ):

    # def count_partitions_combinations(self, var, sizes, cases):
    #     gecode = Solver.lookup("gecode")
    #     model = Model()

    # def count_partitions_combinations(self, var, sizes, cases):
    #     """
    #     Fixes: 
    #     (1) element choices: the way we can choose elements to go in the different sized sub-partitions
    #     (2) size choices: the way a sub_partition of given size can be matched with a partition of given size
    #     """
    #     self.lvl += 1
    #     self.log(f"Sizes {sizes} for cases {cases}")
    #     if len(cases) > 1:
    #         case = cases.pop()
    #         rest_case_sizes = case_sizes.copy()
    #         rest_cases = cases.copy()
    #         count = Solution(0,[])
    #         n = case.size()
    #         k = self.n_vars
    #         count_partitions = self.integer_k_partitions(n,k)
    #         valid_count_partitions = [cp for cp in count_partitions if self.is_feasible_partition(var, sizes, case, cp, case_sizes)]
    #         for cp in valid_count_partitions:
    #             rest_case_sizes[case] = cp
    #             ex_case_ps = self.exchangeable_classes(cp)
    #             ex_sizes = self.exchangeable_classes(sizes)
    #             for s in ex_case_ps:
    #                 ex_case_ps[s] = len(ex_case_ps[s])
    #             for s in ex_sizes:
    #                 ex_sizes[s] = len(ex_sizes[s])
    #             cp_count = self.solve_case_partition(var, cp, case, ex_case_ps, sizes, ex_sizes, rest_cases, rest_case_sizes)
    #             self.log(f"Case {cp} has {cp_count} solutions")
    #             count += cp_count
    #     else:
    #         case = cases[0]
    #         if self.is_feasible_partition(var, sizes, case, sizes):
    #             elem_choices = self.count_partition_choices(sizes, case.size())
    #             n = self.count_partition_choices(sizes, case.size())
    #             count = Solution(n, []).with_choices(elem_choices)
    #             self.log(f"Case {case} has {elem_choices} element choices and {n} solutions (total:{count})")
    #         else:
    #             count = Solution(0,[])
    #             self.log(f"Case {case} is unsat for {sizes}")
    #     self.lvl -= 1
    #     return count

    # def solve_case_partition(self, var, cp, case, ex_case_ps, sizes, ex_sizes, rest_cases, rest_case_sizes):
    #     self.log(f"Consider case-partition {cp} for {case}")
    #     keys = list(ex_sizes.keys())
    #     elem_choices = self.count_partition_choices(cp, case.size())
    #     self.log(f"... has {elem_choices} choices of exchangeable elements")
    #     sub_part_count = Solution(0,[])
    #     # ways of associating the sub_partitions to different sizes
    #     distributions = self.distribute_case_over_sizes(ex_sizes, ex_case_ps)
    #     count = Solution(0,[])
    #     for d in distributions:
    #         rc = rest_cases.copy()
    #         size_choices = self.count_size_choices(ex_sizes, cp, d, keys)
    #         new_sizes = {}
    #         for n_elems in d:
    #             for i, n_sets in enumerate(d[n_elems]):
    #                 # check if enough elements
    #                 if keys[i] >= n_elems:
    #                     new_size = keys[i] - n_elems
    #                     if new_size in new_sizes:
    #                         new_sizes[new_size] += n_sets
    #                     else:
    #                         new_sizes[new_size] = n_sets
    #         rest_sizes = []
    #         for ns in new_sizes:
    #             rest_sizes += [ns] * new_sizes[ns]
    #         rest_sizes += [0] * (len(sizes) - len(rest_sizes))
    #         dist_count = self.count_partitions_combinations(var, rest_sizes, rc, rest_case_sizes)
    #         dist_count = dist_count.with_choices(elem_choices)
    #         dist_count = dist_count.with_choices(size_choices)
    #         sub_part_count += dist_count
    #     return sub_part_count

    # def count_partition_choices(self, distribution, dom_size):
    #     k = distribution[0]
    #     if len(distribution) > 1:
    #         if k <= 1:
    #             return self.count_partition_choices(distribution[1:], dom_size)
    #         else:
    #             choices = math.comb(dom_size,k)
    #             return choices * self.count_partition_choices(distribution[1:], dom_size-k)
    #     else:
    #         if k <= 1:
    #             return 1
    #         else:
    #             choices = math.comb(dom_size,k)
    #             return choices

    # def distribute_case_over_sizes(self, sizes, case):
    #     """
    #     Given the exchangeable classes of valid size partitions, returns the possible ways of matching
    #     a sub-partition of given size to one of the valid partitions' sizes
    #     """
    #     # print("\t", case)
    #     k = len(sizes)
    #     keys = list(sizes.keys())
    #     case_size, n = case.popitem()
    #     # ways of matching n sub-partitions of size case_size across k possible partition sizes
    #     distributions = self.integer_k_partitions(n,k)
    #     options = set()
    #     for d in distributions:
    #         perm = itertools.permutations(d)
    #         options = options.union(set(perm))
    #     # print("\t", case_size, options)
    #     choices = []
    #     for opt in options:
    #         # check that we have enough partitions of that size 
    #         feasible = True
    #         for i in range(0,len(opt)):
    #             s = keys[i]
    #             enough_sets = sizes[s] - opt[i] >= 0
    #             feasible = feasible and enough_sets
    #         # print("\t opt:",opt,feasible)
    #         if feasible:
    #             rest_sizes = sizes.copy()
    #             rest_case = case.copy()
    #             # update: we match opt[i] of size case_size with a partition of size[i]
    #             for i in range(0,len(opt)):
    #                 s = keys[i]
    #                 rest_sizes[s] -= opt[i]
    #             choice = {}
    #             choice[case_size] = opt
    #             if len(case) > 0:
    #                 # print("bch", choice)
    #                 # solve w.r.t the other case_sizes with the remaining available partitions
    #                 rest_choices = self.distribute_case_over_sizes(rest_sizes, rest_case)
    #                 for rest_choice in rest_choices:
    #                     rest_choice.update(choice)
    #                 # print("rch", rest_choices)
    #                 choices += rest_choices
    #             else:
    #                 choices.append(choice)
    #     # print("\t\tchs:", choices)
    #     return choices

    # def count_size_choices(self, ex_sizes, case, distribution, keys):
    #     """
    #     Counts the ways we can distribute some sub-partitions across partitions sizes taking into
    #     account exchangeability of sizes.
    #     """
    #     choices = 1
    #     for size in distribution:
    #         for i, d in enumerate(distribution[size]):
    #             n = ex_sizes[keys[i]]
    #             choices *= math.comb(n,d)
    #     return choices

    def count_satisfied(self, property):
        sat = []
        not_sat = []
        maybe = []
        if isinstance(property, DomainFormula):
            for i, v in enumerate(self.vars):
                pdom = property.domain
                if v.domain in pdom:
                    sat.append(i)
                elif v.domain.disjoint(pdom):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        elif isinstance(property, SizeFormula):
            for i, v in enumerate(self.vars):
                if v.size in property:
                    sat.append(i)
                elif v.size & property == SizeFormula("", P.empty()):
                    not_sat.append(i)
                else:
                    maybe.append(i)
        else: # isinstance(property, SizeFormula)
            for i, v in enumerate(self.vars):
                satisfies = v.satisfies(property)
                if satisfies is None:
                    maybe.append(i)
                elif satisfies:
                    sat.append(i)
                else:
                    not_sat.append(i)
            pass

        return (sat, not_sat, maybe)

    def count_sequence(self, var_list=None):
        count = 1
        vars = var_list if var_list is not None else self.vars
        if not self.alt_type:
            self.log("Counting sequences:")
            for v in vars:
                self.log(v)
                count *= v.domain.size()
                sol = Solution(count, self.histogram())
            self.log("\tDomain product: " + str(count))
        else:
            self.log("Counting permutations:")
            ex_classes = self.exchangeable_classes(vars)
            disjoint_classes = self.disjoint(ex_classes)
            if len(ex_classes) == 1:
                dom = vars[0].domain
                s = dom.size()
                extra = s - len(vars)
                if extra >=0:
                    count = math.factorial(s) // math.factorial(extra)
                    self.log(f"Falling factorial: {count}")
                else: # not enough elements in dom to make a n_vars long permutation
                    self.log(f"{len(vars)} different vars for {s} values!")
                    count = 0
                sol = Solution(count, self.histogram())
            elif disjoint_classes:
                self.log("Different classes but disjoint...")
                for v in vars:
                    self.log(v)
                    count *= v.domain.size()
                sol = Solution(count, self.histogram())
                self.log("\tDomain product: " + str(count))
            else:
                self.log("Splitting injectivity...")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

    def count_subsets(self, var_list = None):
        self.log(f"There are {self.fixed_choices} fixed elements")
        vars = var_list if var_list is not None else self.vars
        y = len(vars) - self.fixed_choices
        dom = vars[self.fixed_choices].domain # first of free vars
        x = dom.size()
        if self.alt_type:
            count = math.comb(x+y-1,y)
            sol = Solution(count, self.histogram())
        else:
            ex_classes = self.exchangeable_classes(vars)
            if len(ex_classes) == 1:
                x = x - self.fixed_choices
                count = math.comb(x,y)
                sol = Solution(count, self.histogram())
            else:
                self.log("Splitting injectivity...")
                scv, rcv = self.split_ex_classes(ex_classes)
                sol = self.split_inj(scv, rcv, [], [])
        return sol

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
        df = dformula
        for case in cases:
            n = cases[case]
            inter = df & case
            exclude = inter.take(n)
            if exclude.domain.size()>0:
                df =  df - exclude
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
                if exclude.domain.size()>0:
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

    def histogram(self):
        ex_classes = self.exchangeable_classes()
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
        parts = []
        for p in self.integer_partitions(n):
            if len(p)<k:
                padded = p + [0]*(k-len(p))
                parts.append(padded)
            elif len(p) == k:
                parts.append(p)
        return parts       

    def is_feasible_split(self, n_split, n_rest, scof, rcof):
        """
        Checks if we ask to observe more properties than available variables
        """
        if scof.values.lower > n_split:
            return False
        if rcof.values.lower > n_rest:
            return False
        return True

    def is_feasible_partition(self, sizes):
        """
        Given a variable LiftedSet, checks if the size of partitions are compatible with the size constraint.
        Similarly, given a case of a counting formula, check if there is some obvious constraint unsatisfiable
        """
        sat = True
        for i, size in enumerate(sizes):
            sat = sat & (size in self.vars[i].size.values)
        return sat

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
            return var & LiftedSet("", property, var.source)
        else: # isinstance(propery, CountingFormula)
            any_size = SizeFormula("", portion.closed(0, portion.inf))
            constr = LiftedSet("", any_size, var, [property]) # var as source?
            return var & constr
    
    def relevant_cases_intersection(self, universe, rest_classes):
        first = rest_classes[0]
        if first.domain == universe:
            return self.relevant_cases_intersection(universe, rest_classes[1:])
        combinations = [[first, first.neg()]]
        for dom in rest_classes[1:]:
            # if we have other subsets we can ignore the universe since there is no element
            # in not(universe)
            if not dom.domain == universe:
                comb = [dom, dom.neg()]
                combinations.append(comb)
        combinations = list(itertools.product(*combinations))
        cases = []
        for c in combinations:
            dom_base = c[0]
            if len(c) > 1:
                for dom in c[1:]:
                    dom_base = dom_base & dom
            if dom_base.domain.size()>0:
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
            all = rest_classes[0]
            none = rest_classes[0].neg()
            res = [none] + [all] 
        else:
            res = self.relevant_cases_intersection(split_df.universe, rest_classes)
        self.log(f"\t{res}")
        return res
       
    # def shatter_exchangeable_partitions(self, var_list=None):
    #     """
    #     All partitions are exchangeable but constants are not: shatter into subproblems where constants are exchangeable.
    #     """
    #     vars = var_list if var_list is not None else self.vars
    #     v = vars[0]
    #     n = v.source.size()
    #     count = Solution(0,[])
    #     size_partitions = self.integer_k_partitions(n, len(vars))
    #     valid_size_partitions = [sp for sp in size_partitions if self.is_feasible_partition(v.size, sp)]
    #     if len(v.constraints) == 0:
    #         return self.count_partitions_by_size()
    #     else:
    #         # get relevant/disjoint set-properties
    #         relevant = list(v.relevant())
    #         cases = self.relevant_cases(v.source, relevant)
    #         # for each valid set of partitions' sizes find the combinations of
    #         # independent partitions that sum up to a valid set size
    #         count = Solution(0,[])     
    #         for vsp in valid_size_partitions:
    #             self.log(f"Considering partition sizes {vsp}")
    #             self.lvl += 1
    #             ex_sizes = self.exchangeable_classes(vsp)
    #             ex_cof_cases = self.get_partition_cof_cases(v, vsp)
    #             # vsp_count = self.count_partitions_combinations(v, vsp, vsp_cases)
    #             self.log(f"...has {vsp_count} solutions.")
    #             count += vsp_count
    #             self.lvl -= 1
    #         return count

    # def get_partition_cof_cases(self, var, valid_sizes):
    #     cof_cases = {} 
    #     for cof in var.constraints:
    #         n = cof.formula.size()
    #         k = len(valid_sizes)
    #         part_cases = self.integer_k_partitions(n, k)
    #         valid_partitions = [pc for pc in part_cases if self.is_feasible_partition(cof.values, pc)]
    #         cof_cases[cof.formula] = valid_partitions
    #     return cof_cases
    
    # def match_sub_partitions(self, case, sizes, case_sizes):
    #     ex_sizes = self.exchangeable_classes(sizes)
    #     keys = list(sizes.keys())
    #     for es in ex_sizes:
    #         es[sizes] = len(es[sizes])
    #     k = len(ex_sizes)
    #     distribution = self.integer_k_partitions(case.size(),k)
    #     valid_distribution = []
    #     for d in distribution:
    #         valid = True
    #         for i in range(1,k):    
    #             s = keys[i]
    #             valid = valid & (distribution[i] >= s*ex_sizes[s])
    #         if valid:
    #             valid_distribution.append(d)
    #     return valid_distribution
                
    def solve(self, log=True):
        self.dolog = log
        self.log(self)
        for c in self.choice_f:
            self.apply_choice(c)
        self.count_f = self.compact_cofs(self.count_f)
        if len(self.count_f) !=0:
            count = self.apply_counts()
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
        split_class_count = self.solve_subproblem(scv, self.type, [], split_class_cofs, self.alt_type)
        if split_class_count != 0:
            self.log("Rest class: ")
            rest_classes_count = self.solve_subproblem(rcv, self.type, [], rest_classes_cofs, self.alt_type)
            count = split_class_count * rest_classes_count
            self.log("==========")
        return count

    def solve_subproblem(self, vars, type, choice_constr, count_constr, alt_type):
        self.log(f"\tSubproblem ({type}):")
        vars = [v.copy() for v in vars]
        subproblem = SharpCSP(vars, type, choice_constr, count_constr, alt_type, self.lvl+1)
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
            count = self.solve_subproblem(self.vars, self.type, [], others, self.alt_type)
        return count.with_choices(n_choices)
        
    def split_ex_classes(self, ex_classes):
        """
        Groups the exchangeable classes as left | right by taking one as left and uniting the others
        """
        it_excls = iter(ex_classes)
        split_class = next(it_excls)
        rest_classes = []
        for i in it_excls:
            rest_classes = rest_classes + ex_classes[i]
        return (ex_classes[split_class], rest_classes)

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
                    self.log(f"Case {n_case} {cases}")
                    split_count_formulas = split_class_cofs.copy()
                    # split_inj_formulas = list(map(lambda i: CountingFormula(cases[i], P.singleton(n_case[i])), range(0,len(n_case))))
                    split_inj_formulas = [CountingFormula(cases[i], P.singleton(n_case[i])) for i in range(0,len(n_case))]
                    split_count_formulas = split_count_formulas + split_inj_formulas
                    self.log(f"Case {n} {split_class} are s.t. {split_count_formulas}:")
                    try:
                        split_count_formulas = self.compact_cofs(split_count_formulas)
                        split_class_count = self.solve_subproblem(split_class_vars, self.type, [], split_count_formulas, self.alt_type)
                        if split_class_count.count != 0: #optimization: do not solve rest if split already unsat
                            for c, hst in split_class_count.histograms:
                                partial_count = Solution(c, hst)
                                filtered = self.filter_domains(hst, rest_classes_vars)
                                rest_count = self.solve_subproblem(filtered, self.type, [], rest_classes_cofs, self.alt_type)
                                self.log(f"Split result = {c} * {rest_count.count} = {partial_count * rest_count}")
                                count += partial_count * rest_count
                    except Unsatisfiable:
                        self.log(f"\t is unsat.")
        self.lvl -= 1
        return count

    def split_partitions(self, split_class_vars, rest_classes_vars, split_class_cofs, rest_classes_cofs):
        self.log(f"Split class :")
        count = Solution(0, [])
        return count

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