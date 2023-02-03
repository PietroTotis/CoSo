import os
import signal
import psutil
import argparse
import time
import random
import signal
import clingo
import time
import portion as P
from multiprocessing import Process, Value
from subprocess import Popen, PIPE, TimeoutExpired
from statistics import mean

from .parser import EmptyException, Parser
from .count import Count
from .configuration import CCounting, CSize
from .util import *

TIMEOUT = 300
random.seed(1234)

ops = [">", "<", "<=", ">=", "!=", "="]
TOOLS = os.path.join(ROOT_DIR, "tools")
TESTS = os.path.join(ROOT_DIR, "tests")
BENCHMARKS = os.path.join(TESTS, "benchmarks")
BENCHMARKS_SYNTH = os.path.join(BENCHMARKS, "synthetic")
RESULTS = os.path.join(TESTS, "results")
ASP_TOOLS = os.path.join(TOOLS, "ASP_tools")
SHARP_SAT = os.path.join(TOOLS, "sharpSAT", "build", "Release", "sharpSAT")
CONJURE = os.path.join(TOOLS, "conjure")


#################
## Miscellanea ##
#################


class Result:
    def __init__(self, solver, solution, time):
        self.solver = solver
        self.solution = solution
        self.time = time

    def __repr__(self):
        return f"{self.solver}: {self.solution} in {self.time}s"

    def __str__(self):
        return f"{self.solver}: {self.solution} in {self.time:.3f}s"


class Context:
    def id(self, x):
        return x

    def seq(self, x, y):
        return [x, y]


def killtree(pid):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.kill()
    parent.kill()


def clean_essence_garbage():
    for file in os.listdir(ROOT_DIR):
        if file.startswith(".MINION"):
            # print(f"Cleaning {file}")
            try:
                os.remove(file)
            except OSError:
                pass
    for file in os.listdir(os.path.join(ROOT_DIR, "src")):
        if file.startswith(".MINION"):
            # print(f"Cleaning {file}")
            try:
                os.remove(file)
            except OSError:
                pass


###############################
## Generation and conversion ##
###############################


def get_n_vars(n):
    vars = []
    for c in range(ord("A"), ord("Z")):
        if len(vars) < n:
            vars.append(chr(c))
    return vars


def get_random_complex_dom(n_domains):
    if n_domains < 2:
        return "dom1"
    else:
        type = random.randint(1, 4)
        if type <= 2:
            n_dom = random.randint(2, n_domains)
            dom = f"dom{n_dom}"
        if type == 2:
            dom = f"¬({dom})"
        if type > 2:
            n_dom_left = random.randint(1, n_domains)
            n_dom_right = random.randint(1, n_domains)
            while n_dom_left == n_dom_right:
                n_dom_right = random.randint(1, n_domains)
            neg_left = random.randint(0, 1)
            neg_right = random.randint(0, 1)
            dom_left = f"dom{n_dom_left}"
            dom_right = f"dom{n_dom_right}"
            if neg_left == 1:
                dom_left = f"¬{dom_left}"
            if neg_right == 1:
                dom_right = f"¬{dom_right}"
            if type == 3:
                dom = f"({dom_left}&{dom_right})"
            else:
                dom = f"({dom_left}+{dom_right})"
        return dom


def generate_problem(
    structure,
    domains_upto,
    universe_size,
    struct_size,
    choice_constraints,
    counting_constraints,
):
    """Generate a CoLa problem

    Args:
        structure (str): structure type, one of: sequence, permutation, subset, multisubset, partition, composition
        domains_upto (int): upper bound to the number of domain to be generated
        u_size (int): upper bounds to the number of elements in each domain
        struct_size (int): [description]
        choice_constraints (bool): [description]
        counting_constraints (bool): [description]

    Returns:
        [type]: [description]
    """
    problem = ""
    problem += f"universe u = {{"
    left = universe_size
    multiset = {}
    elems = []
    for i in range(0, universe_size):
        if left > 0:
            name = f"e{i+1}"
            n_copies = 1 + random.randint(0, left // 2)
            multiset[i] = n_copies
            new_elems = [name] * n_copies
            elems += new_elems
            left -= n_copies
    problem += ", ".join(elems)
    problem += "};\n"
    n_domains = random.randint(2, domains_upto)
    # dom_sizes = {}
    for i in range(1, n_domains):
        base = random.randint(1, len(multiset) - 1)
        upper = random.randint(base + 2, len(multiset) + 1)
        # s = 0
        # for j in range(base,upper):
        #     s += multiset[j]
        # dom_sizes[i] = s
        elems = [f"e{j}" for j in range(base, upper)]
        elems_str = ", ".join(elems)
        problem += f"set dom{i} = {{{elems_str}}};\n"

    # n_domains = random.randint(2, domains_upto)
    # problem += f"set indist uni = [1:{universe_size}];\n"
    # for i in range(1,n_domains):
    #     base = random.randint(1,universe_size-1)
    #     upper = random.randint(base+1,universe_size)
    #     indist = random.randint(0,1) == 1
    #     if indist:
    #         problem+= f"set indist dom{i} = [{base}:{upper}];\n"
    #     else:
    #         problem+= f"set dom{i} = [{base}:{upper}];\n"
    if structure in ["partition", "composition"]:
        problem += f"a in {structure}s(u);\n"
    elif structure == "sequence":
        problem += f"a in [|| u];\n"
    elif structure == "permutation":
        problem += f"a in [| u];\n"
    elif structure == "multisubset":
        problem += f"a in {{|| u}};\n"
    elif structure == "subset":
        problem += f"a in {{| u}};\n"
    else:
        raise Exception(f"Unknown structure {structure}")
    struct_op = random.randint(1, 5)  # avoid large structure with few elements
    problem += f"#a {ops[struct_op]} {struct_size};\n"

    if structure in ["partition", "composition"]:
        if choice_constraints and structure == "composition":
            pos = random.randint(1, struct_size)
            n = random.randint(1, universe_size // 4)
            op = ops[random.randint(0, 5)]
            dom = get_random_complex_dom(n_domains - 1)
            problem += f"#a[{pos}] & {dom} {op} {n} ;\n"
        if counting_constraints:
            n1 = random.randint(1, struct_size // 2)
            op1 = ops[random.randint(0, 5)]
            n2 = random.randint(1, universe_size // 4)
            op2 = ops[random.randint(0, 5)]
            dom = get_random_complex_dom(n_domains - 1)
            problem += f"#{{#{dom} {op2} {n2}}} {op1} {n1} ;\n"
            pass
    else:
        if choice_constraints and structure in {"permutation", "sequence"}:
            n_constr = random.randint(1, struct_size // 2)
            for i in range(1, n_constr):
                pos = random.randint(1, struct_size)
                dom = get_random_complex_dom(n_domains - 1)
                problem += f"a[{pos}] = {dom};\n"
        if counting_constraints:
            n_constr = random.randint(1, struct_size // 2)
            for i in range(1, n_constr):
                n = random.randint(1, struct_size)
                dom = get_random_complex_dom(n_domains - 1)
                op = ops[random.randint(0, 5)]
                problem += f"#{dom}&a {op} {n};\n"
    return problem


def dom2asp(label, domain):
    str = ""
    i = 0
    indist_intervals = domain.elements.find(False)
    for atomic_interval in indist_intervals:
        if atomic_interval != P.empty():
            e = domain.get_label(atomic_interval.lower, atomic_interval.lower)
            l, u = interval_closed(atomic_interval)
            n_copies = u - l + 1
            str += f'{label}_{i}("{e}",{n_copies}).\n'
            str += f"{label}(X) :- {label}_{i}(X, _).\n"
            i += 1
    if len(indist_intervals) == 0:
        dist_intervals = domain.elements.find(True)
        for atomic_interval in dist_intervals:
            l, r = interval_closed(atomic_interval)
            str += f"{label}({l}..{r}).\n"
    else:
        dist_intervals = domain.elements.find(True)
        for atomic_interval in dist_intervals:
            for n in portion.iterate(atomic_interval, step=1):
                e = domain.get_label(n, n)
                str += f'{label}_{i}("{e}", 1).\n'
                str += f"{label}(X) :- {label}_{i}(X, _).\n"
                i += 1
    str += f"universe(X) :- {label}(X).\n"
    return str, i


def problem2asp(problem):
    asp = ""
    n_supports = {}
    for lab, dom in problem.domains.items():
        str, n_support = dom2asp(lab, dom)
        n_supports[lab] = n_support
        asp += str
    if problem.configuration.size is None:
        vals = P.closed(1, problem.universe.size())
        problem.configuration.size = CSize("unconstrained", vals)
    sizes = problem.configuration.size.values
    if problem.configuration.size.values.upper == P.inf:
        ub = problem.universe.size() + 1
        sizes = problem.configuration.size.values.replace(upper=ub)
    lenghts = P.iterate(sizes, step=1)
    sequence = problem.configuration.type == "sequence"
    permutation = problem.configuration.type == "permutation"
    subset = problem.configuration.type == "subset"
    multiset = problem.configuration.type == "multisubset"
    composition = problem.configuration.type == "composition"
    partition = problem.configuration.type == "partition"
    universe_is_set = len(problem.universe.elements.find(False)) == 0
    asp_lengths = []
    for l in lenghts:
        asp_length = ""
        if not composition and not partition:
            vars = get_n_vars(l)
            vars_list = ",".join(vars)
            type = problem.configuration.type
            name = f"{type}_guess_{l}"
            asp_length += f"{name}("
            asp_length += vars_list
            domains = []
            new_props = []
            for i, v in enumerate(vars):
                dom = ""
                for pf in problem.pos_constraints:
                    if pf.pos - 1 == i:
                        if dom != "":
                            raise Exception(
                                "can't correctly translate many pos constraints on same position"
                            )
                        dom_str, _ = dom2asp(f"pf_{i}", pf.formula)
                        new_props.append(dom_str)
                        dom = f"pf_{i}({v})"
                if dom == "":
                    dom = f"universe({v})"
                domains.append(dom)
            asp_length += ") :- " + ", ".join(domains)
            if subset or multiset:
                ineq = "<" if subset and universe_is_set else "<="
                # ineq = "<="
                inequalities = [
                    f"{v}{ineq}{vars[i+1]}"
                    for i, v in enumerate(vars)
                    if i < len(vars) - 1
                ]
            else:
                if permutation and universe_is_set:
                    inequalities = [
                        f"{v1}!={v2}"
                        for i, v1 in enumerate(vars)
                        for j, v2 in enumerate(vars)
                        if i < j
                    ]
                else:
                    inequalities = []
            if len(inequalities) > 0:
                asp_length += ", "
                asp_length += ", ".join(inequalities) + ".\n"
            else:
                asp_length += ".\n"
            asp_length += "\n".join(new_props)
            asp_length += f"1{{{type}_{l}({vars_list}):{name}({vars_list})}}1.\n"
            for k in range(0, l):
                pos_vars = ", ".join(["_" if i != k else "X" for i in range(0, l)])
                asp_length += f"used_{l}(X,{k}) :- {type}_{l}({pos_vars}). \n"
            if not universe_is_set and (permutation or subset):
                for lab in n_supports:
                    for i in range(0, n_supports[lab]):
                        asp_length += f":- {lab}_{i}(S,SN), C = #count{{N:used_{l}(S,N)}}, C>SN.\n"

            for i, cf in enumerate(problem.constraints):
                dlab = f"df_{i}"
                if dlab not in n_supports:
                    dom_str, n = dom2asp(dlab, cf.formula)
                    n_supports[dlab] = n
                    asp += dom_str
                vals = P.closed(0, l) - cf.values
                for n in P.iterate(vals, step=1):
                    asp_length += (
                        f":- C = #count{{N:used_{l}(S,N),df_{i}(S)}}, C={n}.\n"
                    )
        elif composition:
            asp += f"int(0..{problem.universe.size()}).\n"
            for i in range(0, l):
                asp_length += f"part({i}).\n"
            for lab in n_supports:
                for i in range(0, n_supports[lab]):
                    asp_length += f"1{{put(E,N,P): int(N), N<=EN}} 1 :- {lab}_{i}(E, EN), part(P).\n"
                    asp_length += (
                        f":- {lab}_{i}(E,EN), #sum{{N,P:put(E,N,P),part(P)}}!=EN.\n"
                    )
            asp_length += ":- part(P), #count{E,N:put(E,N,P), N>0}==0.\n"
            for i in range(0, l):
                for pf in problem.pos_constraints:
                    if pf.pos - 1 == i:
                        for j, cof in enumerate(pf.formula.ccs):
                            dlab = f"df_{i}_{j}"
                            if dlab not in n_supports:
                                dom_str, n = dom2asp(dlab, cof.formula)
                                n_supports[dlab] = n
                                asp += dom_str
                            vals = P.closed(0, l) - cof.values
                            for n in P.iterate(vals, step=1):
                                asp_length += f":-  C=#sum{{N,E:put(E,N,{i}), {dlab}(E)}}, C={n}.\n"
                        size = pf.formula.size.values
                        if size.lower != 1 or size.upper != l:
                            vals = P.closed(0, l) - size
                            for n in P.iterate(vals, step=1):
                                asp_length += (
                                    f":-  C=#sum{{N,E:put(E,N,{i})}}, C={n}.\n"
                                )

            for i, cf_2 in enumerate(problem.constraints):
                cf_1 = cf_2.formula
                if isinstance(cf_2.formula, CCounting):
                    dlab = f"df_{i}"
                    if dlab not in n_supports:
                        dom_str, n = dom2asp(dlab, cf_1.formula)
                        n_supports[dlab] = n
                        asp += dom_str
                else:
                    dlab = "universe"
                lb = (
                    cf_1.values.lower
                    if cf_1.values.left == P.CLOSED
                    else cf_1.values.lower + 1
                )
                ub = max(l + 1, lb + 1)
                vals = (
                    cf_1.values
                    if cf_1.values.upper != P.inf
                    else cf_1.values.replace(upper=ub)
                )
                count_pred = []
                count_vars = []
                for n in P.iterate(vals, step=1):
                    asp_length += f"cf_{i}_{n}(P,S) :-  S=#sum{{N,E:put(E,N,P), {dlab}(E)}}, part(P), S={n}.\n"
                    count_pred.append(f"C{n}=#count{{P:cf_{i}_{n}(P,{n})}}")
                    count_vars.append(f"C{n}")
                asp_length += f"count_{i}(C) :- "
                asp_length += ", ".join(count_pred) + ", C=" + "+".join(count_vars)
                asp_length += ".\n"
                vals = P.closed(0, l) - cf_2.values
                for n in P.iterate(vals, step=1):
                    asp_length += f":- count_{i}({n}).\n"
        else:
            asp_length = "Partition!"
        asp_lengths.append(asp + asp_length)
    return asp_lengths


def dom2essence(lab, domain):
    dom_str = ""
    indist_intervals = domain.elements.find(False)
    copies = {}
    for atomic_interval in indist_intervals:
        e = domain.get_label(atomic_interval.lower, "e_" + str(atomic_interval.lower))
        e = essence_name(e)
        if atomic_interval != P.empty():
            l, u = interval_closed(atomic_interval)
            n_copies = u - l + 1
            copies[e] = n_copies
    dist_intervals = domain.elements.find(True)
    for atomic_interval in dist_intervals:
        for n in portion.iterate(atomic_interval, step=1):
            e = domain.get_label(n, "e_" + str(n))
            e = essence_name(e)
            copies[e] = 1
    entity_list = ", ".join(copies.keys())
    lab = essence_name(lab)
    if domain.universe == domain or domain.formula == "u":
        dom_str += f"letting universe be new type enum {{ {entity_list} }}\n"
        function_list = ", ".join([f"{e} --> {n}" for e, n in copies.items()])
        dom_str += f"letting f_universe be function({function_list})\n"
    else:
        dom_str += f"letting {lab} be {{ {entity_list} }}\n"
    return dom_str


def range2essence(interval, name, ub):
    if interval.upper == P.inf:
        interval = interval.replace(upper=ub + 1)
    range_vals = ",".join([str(i) for i in P.iterate(interval, step=1)])
    let_str = f"letting {name} be {{ {range_vals} }}\n"
    return let_str


def essence_name(lab):
    if lab in [str(i) for i in range(0, 10)]:
        lab = f"e_{lab}"
    if lab == "true":
        lab = "MYtrue"
    if lab == "false":
        lab = "MYfalse"
    lab = lab.replace("∧", "and")
    lab = lab.replace("¬", "not")
    lab = lab.replace("∨", "or")
    return lab


def problem2essence(problem):
    essence = ""
    added_doms = []
    univ_str = dom2essence("universe", problem.universe)
    essence += univ_str
    for lab, dom in problem.domains.items():
        dom_str = dom2essence(lab, dom)
        if dom_str != univ_str:
            added_doms.append(lab)
            essence += dom_str

    if problem.configuration.size is None:
        ub = problem.universe.size() + 1
        sizes = P.closed(1, problem.universe.size())
        problem.configuration.size = CSize("unconstrained", sizes)
    elif problem.configuration.size.values.upper == P.inf:
        ub = problem.universe.size() + 1
        sizes = problem.configuration.size.values.replace(upper=ub)
    else:
        sizes = problem.configuration.size.values
    lenghts = P.iterate(sizes, step=1)
    sequence = problem.configuration.type == "sequence"
    permutation = problem.configuration.type == "permutation"
    subset = problem.configuration.type == "subset"
    multiset = problem.configuration.type == "multisubset"
    composition = problem.configuration.type == "composition"
    partition = problem.configuration.type == "partition"
    uni = problem.universe.name
    essence_lengths = []
    for l in lenghts:
        name = f"conf_{l}"
        essence_l = f"letting l_{l} be {l}\n"
        constraints = []
        constraint_str = ""
        if not composition and not partition:
            if permutation or subset:
                mset_constraint = f"\tforAll e: {uni}.\n"
                if permutation:
                    mset_constraint += (
                        f"\t\tsum([1 | i: int(1..l_{l}), {name}(i)=e]) <= f_{uni}(e)\n "
                    )
                else:
                    mset_constraint += (
                        f"\t\tforAll e: universe. freq({name},e) <= f_{uni}(e)"
                    )
                constraints.append(mset_constraint)
            for i in range(0, l):
                dom = ""
                for j, pf in enumerate(problem.pos_constraints):
                    if pf.pos - 1 == i:
                        dlab = f"pf_{i}_{j}"
                        if dlab not in added_doms:
                            dom_str = dom2essence(dlab, pf.formula)
                            essence += dom_str
                            added_doms.append(dlab)
                        constraints.append(f"{name}({i+1}) in {dlab}")
            for i, cf in enumerate(problem.constraints):
                dlab = f"df_{i}"
                if dlab not in added_doms:
                    dom_str = dom2essence(dlab, cf.formula)
                    essence += dom_str
                    added_doms.append(dlab)
                if cf.values.upper != P.inf:
                    vals = cf.values
                else:
                    vals = cf.values.replace(upper=l + 1)
                range_vals = ",".join([str(i) for i in P.iterate(vals, step=1)])
                essence_l += f"letting vals_{i} be {{ {range_vals} }}\n"
                if sequence or permutation:
                    constraints.append(
                        f"sum([1 | i: int(1..l_{l}), {name}(i) in {dlab}]) in vals_{i}"
                    )
                else:
                    constraints.append(
                        f"sum([freq({name}, i) | i: {uni}, i in {dlab}]) in vals_{i}"
                    )
            if sequence or permutation:
                essence_l += f"find {name} : sequence (size l_{l}) of {uni}\n"
            if multiset or subset:
                essence_l += f"find {name} : mset (size l_{l}) of {uni}\n"
            if len(constraints) > 0:
                essence_l += "such that \n"
        elif composition:
            myparts = [f"p{i}" for i in range(1, l + 1)]
            essence_l += (
                "letting myparts be new type enum {" + ",".join(myparts) + "}\n"
            )
            essence_l += f"letting n be {problem.universe.size()}\n"
            nonempty = "forAll p: myparts.\n\t sum([put[e,p] | e:universe]) > 0"
            constraints.append(nonempty)
            alldistributed = (
                "forAll e: universe.\n\t sum([put[e,p] | p: myparts]) = f_universe(e)"
            )
            constraints.append(alldistributed)
            ub = problem.universe.size() - l + 1
            for i in range(1, l + 1):
                dom = ""
                for j, pf in enumerate(problem.pos_constraints):
                    if pf.pos == i:
                        if pf.formula.size.values != P.closed(
                            1, ub
                        ) and pf.formula.size.values != P.closed(1, P.inf):
                            name = f"s_{i}"
                            essence_l += range2essence(pf.formula.size.values, name, ub)
                            size_constr = f"sum([put[e,p{i}] | e:{uni}]) in {name}"
                            constraints.append(size_constr)
                        else:
                            for k, cof in enumerate(pf.formula.ccs):
                                df = cof.formula
                                dlab = f"df_{i}_{j}_{k}"
                                if dlab not in added_doms:
                                    dlab = f"df_{i}_{j}_{k}"
                                    dom_str = dom2essence(dlab, df)
                                    essence += dom_str
                                    added_doms.append(dlab)
                                range_name = f"vals_{i}_{j}_{k}"
                                essence_l += range2essence(cof.values, range_name, ub)
                                constraints.append(
                                    f"sum([put[e,p{i}] | e<-{dlab}]) in {range_name}"
                                )
            for i, cf in enumerate(problem.constraints):
                outer_range = range2essence(cf.values, f"vals_{i}_out", l)
                inner_range = range2essence(cf.formula.values, f"vals_{i}_in", ub)
                essence_l += outer_range + inner_range
                dlab = f"df_{i}"
                if isinstance(cf.formula, CSize):
                    dom_str = "universe"
                else:
                    if dlab not in added_doms:
                        dom_str = dom2essence(dlab, cf.formula.formula)
                        essence += dom_str
                        added_doms.append(dlab)
                cf_constraint = f"|[p | p:myparts, sum([put[e,p] | e <- {dom_str}]) in vals_{i}_in]| in vals_{i}_out"
                constraints.append(cf_constraint)

            essence_l += (
                f"find put: matrix indexed by [universe, myparts] of int(0..n)\n"
            )
            essence_l += "such that \n"
        else:
            essence_l += ""
        constraint_str = "\n\t/\ ".join(constraints)
        essence_l += constraint_str
        essence_l += "\n"
        essence_l = essence + essence_l
        essence_lengths.append(essence_l)
    return essence_lengths


def problem2cnf(problem):
    programs = problem2asp(problem)
    gringo = os.path.join(ASP_TOOLS, "gringo")
    lp2normal = os.path.join(ASP_TOOLS, "lp2normal-2.18")
    lp2sat = os.path.join(ASP_TOOLS, "lp2sat-1.24")
    lp2atomic = os.path.join(ASP_TOOLS, "lp2atomic-1.17")
    input = os.path.join(ASP_TOOLS, "tmp.lp")
    out = os.path.join(ASP_TOOLS, "out.cnf")
    gringo_out = os.path.join(ASP_TOOLS, "gringo_out.lp")
    nor_out = os.path.join(ASP_TOOLS, "nor_out.lp")
    at_out = os.path.join(ASP_TOOLS, "at_out.cnf")

    cnf_lengths = []

    for program in programs:
        lp = open(input, "w+")
        lp.write(program)
        lp.close()
        try:
            p = Popen(
                [f"{gringo} {input} > {gringo_out}"],
                shell=True,
            )
            p.wait()
            p = Popen([f"{lp2normal} {gringo_out} > {nor_out}"], shell=True)
            p.wait()
            p = Popen([f"{lp2atomic} {nor_out} > {at_out}"], shell=True)
            p.wait()
            p = Popen([f"{lp2sat} {at_out} > {out}"], shell=True)
            p.wait()

            translation = open(out, "r").read()
            cnf_lengths.append(translation)

        except Exception:
            p.terminate()
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        # print(n_program)
    os.remove(input)
    os.remove(out)
    os.remove(gringo_out)
    os.remove(nor_out)
    os.remove(at_out)
    return cnf_lengths


def generate_constrained(folder, pconstr, cconstr):
    size_domain = [10, 15, 20]
    size_struct = [5, 10, 15]
    cases = [
        "sequence",
        "permutation",
        "subset",
        "multisubset",
        "partition",
        "composition",
    ]
    for d in size_domain:
        domains_upto = d // 3
        for s in size_struct:
            for case in cases:
                case_path = os.path.join(folder, case)
                os.makedirs(case_path, exist_ok=True)
                if d > s:
                    # print(f"{name}, size domain={d}, size structure={s}")
                    problem = generate_problem(
                        case, domains_upto, d, s, pconstr, cconstr
                    )
                    filename = f"{case}_{d}_{s}.test"
                    path = os.path.join(case_path, filename)
                    file = open(path, "w")
                    file.write(problem)
                    file.close()
                    parser = Parser(str(path))
                    parser.parse()
                    # problem_asp = problem2asp(parser.problem)
                    # filename_asp = f"{case}_{d}_{s}.lp"
                    # path_asp = os.path.join(case_path, filename_asp)
                    # file_asp = open(path_asp, "w")
                    # file_asp.write(problem_asp)
                    # file_asp.close()


#####################
## Running solvers ##
#####################


def run_asp(programs, count):
    n = 0
    for program in programs:
        # print(program)
        ctl = clingo.Control()
        ctl.configuration.solve.models = 0
        ctl.add("base", [], program)
        ctl.ground([("base", [])], context=Context())
        with ctl.solve(yield_=True, async_=True) as handle:
            n_length = sum(1 for _ in handle)  # +1
            n += n_length
        count.value = n


def run_sat(programs, count):
    count.value = 0
    gringo = os.path.join(ASP_TOOLS, "gringo")
    lp2normal = os.path.join(ASP_TOOLS, "lp2normal-2.18")
    lp2sat = os.path.join(ASP_TOOLS, "lp2sat-1.24")
    lp2atomic = os.path.join(ASP_TOOLS, "lp2atomic-1.17")
    input = os.path.join(ASP_TOOLS, "tmp.lp")
    out = os.path.join(ASP_TOOLS, "out.cnf")
    gringo_out = os.path.join(ASP_TOOLS, "gringo_out.lp")
    nor_out = os.path.join(ASP_TOOLS, "nor_out.lp")
    at_out = os.path.join(ASP_TOOLS, "at_out.cnf")
    for program in programs:
        lp = open(input, "w+")
        lp.write(program)
        lp.close()
        try:
            p = Popen(
                [f"{gringo} {input} > {gringo_out}"],
                shell=True,
            )
            p.wait()
            p = Popen([f"{lp2normal} {gringo_out} > {nor_out}"], shell=True)
            p.wait()
            p = Popen([f"{lp2atomic} {nor_out} > {at_out}"], shell=True)
            p.wait()
            p = Popen([f"{lp2sat} {at_out} > {out}"], shell=True)
            p.wait()

            p = Popen(
                [f"{SHARP_SAT} -t {TIMEOUT} {out}"],
                shell=True,
                stdout=PIPE,
                stderr=PIPE,
                start_new_session=True,
            )
            std_out, std_err = p.communicate(timeout=TIMEOUT)
            sol = std_out.decode("UTF-8")
            n_string = (
                sol[sol.find("# solutions") + 11 : sol.find("# END")]
                .replace("\n", "")
                .replace(" ", "")
            )
            n_program = int(n_string)
            count.value += n_program
        except Exception:
            p.terminate()
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            count.value = -1
            break
        # print(n_program)
    os.remove(input)
    os.remove(out)
    os.remove(gringo_out)
    os.remove(nor_out)
    os.remove(at_out)


def call_essence(exec_conjure, conjure_output, input, conjure_env):
    p = Popen(
        [
            # f"{exec_conjure} solve -ac -o {conjure_output} --solutions-in-one-file --number-of-solutions=all --limit-time {TIMEOUT}  --log-level lognone {input}"
            f"{exec_conjure} solve -ac -o {conjure_output} --solutions-in-one-file --number-of-solutions=all --log-level lognone {input}"
        ],
        start_new_session=True,
        shell=True,
        env=conjure_env,
    )
    p.wait()


def run_essence(programs, count):
    count.value = 0
    for program in programs:
        exec_conjure = os.path.join(CONJURE, "conjure")
        input = os.path.join(TOOLS, "model.essence")
        conjure_output = os.path.join(TOOLS, "conjure-output")
        output = os.path.join(conjure_output, "model000001.solutions")
        model = open(input, "w+")
        model.write(program)
        model.close()
        conjure_env = os.environ.copy()
        conjure_env["PATH"] = os.path.abspath(CONJURE) + ":" + conjure_env["PATH"]
        try:
            p = Process(
                target=call_essence,
                args=(exec_conjure, conjure_output, input, conjure_env),
            )
            p.start()
            p.join(TIMEOUT)
            n = 0
            if p.is_alive():
                killtree(p.pid)
                p.join()
                print("Killed")
            else:
                if os.path.exists(output):
                    with open(output, "r") as sol:
                        n_prog = sol.read().count("Count:")
                        n += n_prog
            count.value += n
        except Exception:
            p.terminate()
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            count.value = -1
            break
        finally:
            clean_essence_garbage()
    clean_essence_garbage()


def run_coso_timeout(problem, count, n_subproblems, log=False):
    sol = problem.solve(log)
    count.value = sol.count
    n_subproblems.value = sol.subproblems


def run_coso(problem, log=False):
    print("Running CoSo...")
    try:
        count = Value("i", 0, lock=False)
        n_subproblems = Value("i", 0, lock=False)
        start = time.time()
        p = Process(target=run_coso_timeout, args=(problem, count, n_subproblems, log))
        start = time.perf_counter()
        p.start()
        p.join(TIMEOUT)
        finish = time.perf_counter()
        if p.is_alive():
            killtree(p.pid)
            p.join()
            print("Killed")
            res = Result("CoSo", Count(-1, [], None, -1), TIMEOUT)
        else:
            sol = Count(count.value, [], None, n_subproblems.value)
            res = Result("CoSo", sol, finish - start)
            print(res)
    except:
        res = Result("CoSo", Count(-1, [], None, -1), TIMEOUT)
        print(f"CoSo timeout")
    return res


def run_solver(problem, sys_name, translate, run):
    t = translate(problem)
    print(f"Running {sys_name}...")
    try:
        count = Value("i", 0, lock=False)
        start = time.time()
        p = Process(target=run, args=(t, count))
        start = time.perf_counter()
        p.start()
        p.join(TIMEOUT)
        finish = time.perf_counter()
        if p.is_alive():
            killtree(p.pid)
            p.join()
            print("Killed")
            res = Result(sys_name, -1, TIMEOUT)
        elif p.exitcode != 0:
            res = Result(sys_name, -1, TIMEOUT)
        else:
            res = Result(sys_name, count.value, finish - start)
            print(res)
    except:
        killtree(p.pid)
        res = Result(sys_name, -1, TIMEOUT)
        print(f"{sys_name} timeout")
    finally:
        if sys_name == "Essence":
            clean_essence_garbage()
    return res


def compare(problem, name, sys_name, translate, run, log=False):
    print(f"Comparing on {name}")
    # run_coso(problem, log)
    run_solver(problem, sys_name, translate, run)


def count_sols(filename):
    f = open(filename, "r")
    n = 0
    for l in f.read():
        if l[-1] == ";":
            n += 1
    f.close()
    return n


##################
## Benchmarking ##
##################


def test_folder(folder, asp, sat, essence, start_from=None):
    results_coso = {}
    results_asp = {}
    results_sat = {}
    results_essence = {}
    if start_from is not None:
        run = False
    else:
        run = True
    for filename in os.listdir(folder):
        run = run or (filename == start_from)
        if run and (filename.endswith(".test") or filename.endswith(".pl")):
            print(f"Test {filename}:")
            parser = Parser(os.path.join(folder, filename))
            try:
                parser.parse()
                problem = parser.problem
                partition = problem.configuration.type == "partition"
                if not partition:
                    if essence:
                        res = run_solver(
                            parser.problem, "Essence", problem2essence, run_essence
                        )
                        results_essence[filename] = res
                    if asp:
                        res = run_solver(parser.problem, "Clingo", problem2asp, run_asp)
                        results_asp[filename] = res
                    if sat:
                        res = run_solver(
                            parser.problem, "SharpSAT", problem2asp, run_sat
                        )
                        results_sat[filename] = res
                res = run_coso(problem)
                results_coso[filename] = res
            except EmptyException:
                print("Empty: skipping")

    folder_name = os.path.basename(folder)
    results_file = f"bench_results_{folder_name}.csv"
    with open(os.path.join(RESULTS, results_file), "w+") as f:
        f.write(f"benchmark\tsolver\tn_subproblems\tn_solutions\ttime\n")
        f.close()
    if len(results_coso) > 0:
        export_results(results_coso, results_file, "CoSo")
        with open(os.path.join(RESULTS, f"subs_results_{folder_name}.csv"), "w+") as f:
            export_coso_results(results_coso, f)
    if len(results_asp) > 0:
        export_results(results_asp, results_file, "ASP")
    if len(results_sat) > 0:
        export_results(results_sat, results_file, "#SAT")
    if len(results_essence) > 0:
        export_results(results_essence, results_file, "Essence")


def export_results(results, file, solver):
    present_results(results, solver)
    with open(os.path.join(RESULTS, file), "a") as f:
        for res in results:
            if solver == "CoSo":
                n_sub = results[res].solution.subproblems
                sol = results[res].solution.count
            else:
                n_sub = -1
                sol = results[res].solution
            f.write(f"{res}\t{solver}\t{n_sub}\t{sol}\t{results[res].time}\n")
    f.close()


def coso_subproblem_stats(results):
    groupby_subs_time = {}
    groupby_subs_num = {}
    for name in results:
        n = results[name].solution.subproblems
        t = results[name].time
        if n in groupby_subs_time:
            groupby_subs_time[n] = mean([groupby_subs_time[n], t])
            groupby_subs_num[n] += 1
        else:
            groupby_subs_num[n] = 1
            groupby_subs_time[n] = t
    return (groupby_subs_num, groupby_subs_time)


def present_results(results, solver):
    times = [results[name].time for name in results]
    avg_time = mean(times)
    print(f"Average {solver} time: {avg_time:.2f}")
    if solver == "CoSo":
        gsn, gst = coso_subproblem_stats(results)
        for n in gsn:
            print(f"{n} subs (n={gsn[n]}): {gst[n]:.2f}s")


def export_coso_results(results, file):
    gsn, gst = coso_subproblem_stats(results)
    file.write("n_subproblems\tcount\ttime\n")
    for n in gsn:
        file.write(f"{n}\t{gsn[n]}\t{gst[n]}\n")
    file.close()


def run_benchmarks(plot, start_from):
    types = []
    # types = ["subset", "composition", "sequence"]
    # types = ["multisubset", "permutation", "sequence", "subset", "composition"]
    for type in types:
        dir = os.path.join(BENCHMARKS_SYNTH, type)
        print(dir)
        test_folder(dir, True, True, True, start_from)
        clean_essence_garbage()

    examples = os.path.join(TESTS, "examples")
    test_folder(examples, True, False, False, start_from)

    # grow_doms = os.path.join(BENCHMARKS, "growing_domains")
    # results_asp = {}
    # results_sat = {}
    # results_essence = {}
    # test_folder(grow_doms, False, False, False)
    # asp_dir = os.path.join(grow_doms, "ASP")
    # essence_dir = os.path.join(grow_doms, "Essence")
    # translate = lambda x: [str(x)]
    # for asp_problem in os.listdir(asp_dir):
    #     with open(os.path.join(asp_dir, asp_problem), "r") as f:
    #         problem = f.read()
    #         f.close()
    #         res_asp = run_solver(problem, "Clingo", translate, run_asp)
    #         res_sat = run_solver(problem, "SharpSAT", translate, run_sat)
    #         results_asp[asp_problem] = res_asp
    #         results_sat[asp_problem] = res_sat

    # for essence_problem in sorted(os.listdir(os.path.join(grow_doms, "Essence"))):
    #     with open(os.path.join(essence_dir, essence_problem), "r") as f:
    #         problem = f.read()
    #         f.close()
    #         res_essence = run_solver(problem, "Conjure", translate, run_essence)
    #         results_essence[essence_problem] = res_essence
    # results_file = f"bench_results_growing_domains.csv"
    # clean_essence_garbage()
    # if len(results_asp) > 0:
    #     export_results(results_asp, results_file, "ASP")
    # if len(results_sat) > 0:
    #     export_results(results_sat, results_file, "#SAT")
    # if len(results_essence) > 0:
    #     export_results(results_essence, results_file, "Essence")
    # plot(False)


def plot(export):
    from gen_plots import plot as plot_results

    plot_results(BENCHMARKS_SYNTH, pgf=export)


def translate_folder(folder, translation, extension=None):
    for filename in os.listdir(folder):
        if filename.endswith(".test") or filename.endswith(".pl"):
            path = os.path.join(folder, filename)
            print(f"Translating {filename}:")
            translate(path, translation, extension)


def translate(file, translation, tool=None):
    name = os.path.basename(file).split(".")[0]
    if tool == "Clingo":
        extension = "lp"
    elif tool == "SharpSAT":
        extension = "cnf"
    elif tool == "Essence":
        extension = "essence"
    else:
        extension = "txt"
    loc = os.path.dirname(file)
    parser = Parser(file)
    parser.parse()
    problem = parser.problem
    problems_translated = translation(problem)
    if len(problems_translated) > 1:
        out = os.path.join(loc, f"{name}.{extension}")
        os.makedirs(out, exist_ok=True)
        for i, pt in enumerate(problems_translated):
            pt_name = name + f"_{i}.txt"
            pt_path = os.path.join(out, pt_name)
            with open(pt_path, "w+") as f:
                f.write(pt)
                f.close()
    else:
        p = problems_translated[0]
        out = os.path.join(loc, f"{name}.{extension}")
        with open(out, "w+") as f:
            f.write(p)
            f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Run solver on file")
    parser.add_argument("-d", help="Run solver on directory")
    parser.add_argument("-b", action="store_true", help="Run benchmarks")
    parser.add_argument("-l", action="store_true", help="Activate logging")
    parser.add_argument("-g", help="Geneate random tests in given folder")
    parser.add_argument(
        "-t", "--timeout", type=int, default=TIMEOUT, help="Set timeout in seconds"
    )
    parser.add_argument(
        "-e",
        action="store_true",
        help="Export plots or translation only without running solvers",
    )
    parser.add_argument("--plot", action="store_true", help="Plot benchmarks results")
    parser.add_argument("--asp", action="store_true", help="Compare with 'asp'")
    parser.add_argument("--sat", action="store_true", help="Compare with 'sharpSAT'")
    parser.add_argument("--essence", action="store_true", help="Compare with 'essence'")
    parser.add_argument(
        "--sf", help="Start running benchmark in folder after given filename"
    )
    parser.add_argument(
        "--noposconstr",
        action="store_false",
        help="Disable positional constraint generation when -g",
    )
    parser.add_argument(
        "--nocountconstr",
        action="store_false",
        help="Disable counting constraint generation when -g",
    )
    args = parser.parse_args()
    TIMEOUT = args.timeout

    solvers = []
    if args.asp:
        solvers.append(("Clingo", problem2asp, run_asp))
    if args.sat:
        solvers.append(("SharpSAT", problem2cnf, run_sat))
    if args.essence:
        solvers.append(("Essence", problem2essence, run_essence))
    if args.f:
        if args.e:
            for sname, translation, _ in solvers:
                translate(args.f, translation, sname)
        else:
            parser = Parser(args.f)
            parser.parse()
            run_coso(parser.problem, args.l)
            for solver_args in solvers:
                compare(parser.problem, args.f, *solver_args, log=args.l)
    elif args.d:
        if args.e:
            for sname, translation, _ in solvers:
                translate_folder(args.d, translation, sname)
        else:
            test_folder(args.d, args.asp, args.sat, args.essence)
    elif args.b:
        run_benchmarks(args.plot, args.sf)
    elif args.g:
        generate_constrained(args.g, args.noposconstr, args.nocountconstr)
    elif args.plot:
        plot(args.e)
    else:
        parser.print_help()
