import errno
import os
import signal
import functools
import argparse
import time
import portion
import pyinotify
import random
import signal
import clingo
import portion as P
from subprocess import Popen, PIPE
from parser import EmptyException, Parser
from formulas import PosFormula, And, Or, Not
from util import interval_closed, ROOT_DIR

TIMEOUT = 20
random.seed(123)

class Context:
    def id(self, x):
        return x
    def seq(self, x, y):
        return [x, y]


class TimeoutError(Exception):
    pass

def timeout(seconds=TIMEOUT, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator

ops = [">","<","<=",">=","!=","="]
tools = os.path.join(ROOT_DIR, "tools")
asp_tools = os.path.join(tools, "ASP_tools")
sharp_sat = os.path.join(tools, "sharpSAT", "build", "Release", "sharpSAT")
conjure = os.path.join(tools, "conjure")

class ModHandler(pyinotify.ProcessEvent):
    def process_IN_CLOSE_WRITE(self, evt):
        pass


def get_n_vars(n):
    vars = []
    for c in range(ord('A'),ord('Z')):
        if len(vars) < n:
            vars.append(chr(c))
    return vars

def get_random_complex_dom(n_domains):
    if n_domains < 2:
        return "dom1"
    else:
        type = random.randint(1,4)
        if type <= 2:
            n_dom = random.randint(2,n_domains)
            dom = f"dom{n_dom}"
        if type == 2:
            dom = f"¬({dom})"
        if type > 2:
            n_dom_left = random.randint(1,n_domains)
            n_dom_right = random.randint(1,n_domains)
            while n_dom_left == n_dom_right:
                n_dom_right = random.randint(1,n_domains)
            neg_left = random.randint(0,1)
            neg_right = random.randint(0,1)
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

def generate_problem(structure,domains_upto,universe_size,struct_size,choice_constraints,counting_constraints):
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
    problem += f"set universe = {{"
    left = universe_size
    multiset = {}
    elems = []
    for i in range(0, universe_size):
        if left > 0:
            name  = f"e{i+1}"
            n_copies = 1+ random.randint(0,left//2)
            multiset[i] = n_copies
            new_elems = [name] * n_copies
            elems += new_elems
            left -= n_copies
    problem += ", ".join(elems)
    problem += "};\n"           
    n_domains = random.randint(2, domains_upto)
    # dom_sizes = {}
    for i in range(1,n_domains):
        base = random.randint(1,len(multiset)-1)
        upper = random.randint(base+2,len(multiset)+1)
        # s = 0
        # for j in range(base,upper):
        #     s += multiset[j]
        # dom_sizes[i] = s
        elems = [f"e{j}" for j in range(base,upper)]
        elems_str = ", ".join(elems)
        problem+= f"set dom{i} = {{{elems_str}}};\n"

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
        problem += f"a in {structure}s(universe);\n"
    elif structure == "sequence":
        problem += f"a in [|| universe];\n"
    elif structure == "permutation":
        problem += f"a in [| universe];\n"
    elif structure == "multisubset":
        problem += f"a in {{|| universe}};\n"
    elif structure == "subset":
        problem += f"a in {{| universe}};\n"
    else:
        raise Exception(f"Unknown structure {structure}")
    struct_op = random.randint(1,5) # avoid large structure with few elements
    problem += f"#a {ops[struct_op]} {struct_size};\n"

    if structure in ["partition", "composition"]:
        if choice_constraints and structure == "composition":
            pos = random.randint(1,struct_size)
            n = random.randint(1,universe_size // 4)
            op = ops[random.randint(0,5)]
            dom = get_random_complex_dom(n_domains-1)
            problem += f"#a[{pos}] & {dom} {op} {n} ;\n"
        if counting_constraints:
            n1 = random.randint(1,struct_size // 2)
            op1 = ops[random.randint(0,5)]
            n2 = random.randint(1,universe_size // 4 )
            op2 = ops[random.randint(0,5)]
            dom = get_random_complex_dom(n_domains-1)
            problem += f"#{{#{dom} {op2} {n2}}} {op1} {n1} ;\n"
            pass
    else:
        if choice_constraints and structure in {"permutation", "sequence"}:
            n_constr = random.randint(1,struct_size // 2)
            for i in range(1,n_constr):
                pos = random.randint(1,struct_size)
                dom = get_random_complex_dom(n_domains-1)
                problem += f"a[{pos}] = {dom};\n"
        if counting_constraints:
            n_constr = random.randint(1,struct_size // 2)
            for i in range(1,n_constr):
                n = random.randint(1,struct_size)
                dom = get_random_complex_dom(n_domains-1)
                op = ops[random.randint(0,5)]
                problem += f"#{dom}&a {op} {n};\n"
    return problem

def domf2minizinc(df,n):
    if isinstance(df,Not):
        aux_d, name, n = domf2minizinc(df.child,n)
        descr = f"\nset of int: df{n+1}= uni diff {name};\n"
        return aux_d + descr, f"df{n+1}", n+1
    elif isinstance(df,And):
        aux_dl, namel, n = domf2minizinc(df.left,n)
        aux_dr, namer, m = domf2minizinc(df.right,n)
        descr = f"\nset of int: df{m+1}= {namel} intersect {namer};\n"
        return aux_dl+ aux_dr + descr, f"df{m+1}", m+1
    elif isinstance(df, Or):
        aux_dl, namel, n = domf2minizinc(df.left,n)
        aux_dr, namer, m = domf2minizinc(df.right,n)
        descr = f"\nset of int: df{m+1}= {namel} union {namer};\n"
        return aux_dl+ aux_dr + descr, f"df{m+1}", m+1
    else:
        return  "", df, n

def problem2minizinc(problem):
    minizinc = 'include "globals.mzn";\n'
    for dom in problem.domains.values():
        minizinc += f"set of int: {dom.name} = {{"
        minizinc += ",".join(list(map(str,portion.iterate(dom.elements.domain(), step =1))))
        minizinc += "};\n"
    length = problem.structure.size.values.lower
    minizinc += f"int: n = {length};\n"
    sequence = problem.structure.type == "sequence" or problem.structure.type == "permutation"
    subset = problem.structure.type == "subset" or problem.structure.type == "multisubset"
    alldiff = problem.structure.type == "permutation" or problem.structure.type == "subset"
    if sequence:
        minizinc += f"array[1..n] of var uni: sequence;\n"
        if alldiff:
            minizinc += "constraint alldifferent(sequence);\n"
    elif subset:
        minizinc += f"var set of uni: sub;\n"
        minizinc += f"constraint card(sub) == n;\n"
    n = 0
    for chf in problem.choice_formulas:
        if isinstance(chf, PosFormula):
            aux_doms, dom_f, n  = domf2minizinc(chf.dformula.name, n)
            minizinc += aux_doms
            minizinc += f"constraint sequence[{chf.pos}] in {dom_f};\n"
        elif isinstance(chf, InFormula):
            minizinc += f"constraint {chf.entity} in sets;\n"
    for cof in problem.count_formulas:
        if sequence:
            range = "i in 1..n"
            elem = "sequence[i]"
        elif subset:
            range = "i in sub"
            elem = "i"
        intv = cof.values
        bounds = []
        for left, lower, upper, right in portion.to_data(intv):
            if upper > 1000:#some issues comparing to portion.inf
                upper = length
            if not left:
                lower+= 1
            if not right:
                upper+= 1
            bounds.append((lower,upper))
        aux_doms, dom_f, n  = domf2minizinc(cof.formula.name, n)
        for l in aux_doms.split('\n'):
            if l not in minizinc:
                minizinc += l + "\n"
        minizinc += f"constraint "
        cofs = []
        for lower, upper in bounds:
            cofs.append(f"sum({range})(bool2int({elem} in {dom_f})) >= {lower} /\\ sum({range})(bool2int({elem} in {dom_f})) <= {upper}")
        minizinc += "\\/\n".join(cofs) + ";\n"
    minizinc += "solve satisfy;"
    return minizinc

def dom2asp(label, domain):
    str = ""
    i = 0
    indist_intervals = domain.elements.find(False)
    for atomic_interval in indist_intervals:
        e = domain.labels.get(atomic_interval.lower,atomic_interval.lower)
        if atomic_interval != P.empty():
            l, u = interval_closed(atomic_interval)
            n_copies = u-l+1
            str += f"{label}_{i}(\"{e}\",{n_copies}).\n"
            str += f"{label}(X) :- {label}_{i}(X, _).\n"
            i += 1
    dist_intervals = domain.elements.find(True)
    for atomic_interval in dist_intervals:
        for e in portion.iterate(atomic_interval, step =1):
            e = domain.labels.get(atomic_interval.lower,atomic_interval.lower)
            str += f"{label}_{i}(\"{e}\", 1).\n"
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
    sizes = problem.configuration.size.values
    if problem.configuration.size.values.upper == P.inf:
        ub = problem.universe.size() +1
        sizes = problem.configuration.size.values.replace(upper = ub)
    lenghts = P.iterate(sizes, step=1)
    sequence = problem.configuration.type == "sequence" 
    permutation =  problem.configuration.type == "permutation"
    subset = problem.configuration.type == "subset" 
    multiset =  problem.configuration.type == "multisubset"
    composition =  problem.configuration.type == "composition"
    partition =  problem.configuration.type == "partition"
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
                for pf in problem.pos_formulas:
                    if pf.pos-1 == i:
                        if dom != "":
                            raise Exception("can't correctly translate many pos constraints on same position")
                        dom_str, _ = dom2asp(f"pf_{i}", pf.formula)
                        new_props.append(dom_str)
                        dom = f"pf_{i}({v})"
                if dom == "":
                    dom = f"universe({v})" 
                domains.append(dom)
            asp_length += ") :- " + ", ".join(domains)
            if subset or multiset:
                # ineq = "<" if subset else "<="
                ineq = "<="
                inequalities = [f"{v}{ineq}{vars[i+1]}" for i, v in enumerate(vars) if i < len(vars)-1]
                if len(inequalities) > 0:
                    asp_length += ", "
                    asp_length += ", ".join(inequalities) + ".\n"
                else:
                    asp_length += ".\n"
            else:
                asp_length += ".\n"
            asp_length += "\n".join(new_props)
            asp_length += f"1{{{type}_{l}({vars_list}):{name}({vars_list})}}1.\n"
            for k in range(0,l):
                pos_vars = ", ".join(["_" if i != k else "X" for i in range(0,l) ])
                asp_length += f"used_{l}(X,{k}) :- {type}_{l}({pos_vars}). \n"
            if permutation or subset:
                for lab in n_supports:
                    for i in range(0,n_supports[lab]):
                        asp_length += f":- {lab}_{i}(S,SN), C = #count{{N:used_{l}(S,N)}}, C>SN.\n"
            
            for i, cf in enumerate(problem.count_formulas):
                dlab = f"df_{i}"
                if dlab not in n_supports:
                    dom_str, n = dom2asp(dlab, cf.formula)
                    n_supports[dlab] = n
                    asp += dom_str
                vals = P.closed(0,l) - cf.values
                for n in P.iterate(vals, step=1):
                    asp_length += f":- C = #count{{N:used_{l}(S,N),df_{i}(S)}}, C={n}.\n"
        elif composition:
            asp += f"int(0..{problem.universe.size()}).\n"
            for i in range(0,l):
                asp_length += f"part({i}).\n" 
            for lab in n_supports:
                for i in range(0,n_supports[lab]):
                    asp_length += f"1{{put(E,N,P): int(N), N<=EN}} 1 :- {lab}_{i}(E, EN), part(P).\n"
                    asp_length += f":- {lab}_{i}(E,EN), #sum{{N,P:put(E,N,P),part(P)}}!=EN.\n"
            asp_length += ":- part(P), #count{E,N:put(E,N,P), N>0}==0.\n"
            for i in range(0,l):
                for pf in problem.pos_formulas:
                    if pf.pos-1 == i:
                        for j, cof in enumerate(pf.formula.cofs):
                            dlab = f"df_{i}_{j}"
                            if dlab not in n_supports:
                                dom_str, n = dom2asp(dlab, cof.formula)
                                n_supports[dlab] = n
                                asp += dom_str
                            vals = P.closed(0,l) - cof.values
                            for n in P.iterate(vals, step=1):
                                asp_length += f":-  C=#sum{{N,E:put(E,N,{i}), {dlab}(E)}}, C={n}.\n"
                        size = pf.formula.size.values
                        if size.lower != 1 or size.upper != l:
                            vals = P.closed(0,l) - size
                            for n in P.iterate(vals, step=1):
                                asp_length += f":-  C=#sum{{N,E:put(E,N,{i})}}, C={n}.\n"
                        
            for i, cf_2 in enumerate(problem.count_formulas):
                cf_1 = cf_2.formula
                dlab = f"df_{i}"
                if dlab not in n_supports:
                    dom_str, n = dom2asp(dlab, cf_1.formula)
                    n_supports[dlab] = n
                    asp += dom_str
                vals = cf_1.values if cf_1.values.upper != P.inf else cf_1.values.replace(upper=l+1)
                count_pred = []
                count_vars = []
                for n in P.iterate(vals, step=1):
                    asp_length += f"cf_{i}_{n}(P,S) :-  S=#sum{{N,E:put(E,N,P), {dlab}(E)}}, part(P), S={n}.\n"
                    count_pred.append(f"C{n}=#count{{P:cf_{i}_{n}(P,{n})}}")
                    count_vars.append(f"C{n}")
                asp_length += f"count_{i}(C) :- "
                asp_length += ", ".join(count_pred) + ", C=" + "+".join(count_vars)
                asp_length += ".\n"
                vals = P.closed(0,l) - cf_2.values
                for n in P.iterate(vals, step=1):
                    asp_length += f":- count_{i}({n}).\n"
        else:
            pass
        asp_lengths.append(asp+asp_length)
    return asp_lengths

def dom2essence(lab, domain):
    dom_str = ""
    indist_intervals = domain.elements.find(False)
    copies = {}
    for atomic_interval in indist_intervals:
        e = domain.labels.get(atomic_interval.lower,atomic_interval.lower)
        if atomic_interval != P.empty():
            l, u = interval_closed(atomic_interval)
            n_copies = u-l+1
            copies[e] = n_copies
    dist_intervals = domain.elements.find(True)
    for atomic_interval in dist_intervals:
        for n in portion.iterate(atomic_interval, step =1):
            e = domain.labels.get(n,n)
            copies[e] = 1
    entity_list = ", ".join(copies.keys())
    if domain.universe == domain:
        dom_str += f"letting {lab} be new type enum {{ {entity_list} }}\n"
    else:
        dom_str += f"letting {lab} be {{ {entity_list} }}\n"
    function_list = ", ".join([f"{e} --> {n}" for e,n in copies.items()])
    dom_str += f"letting f_{lab} be function({function_list})\n"
    return dom_str

def range2essence(interval, name, ub):
    if interval.upper == P.inf:
        interval = interval.replace(upper = ub+1)
    range_vals = ",".join([str(i) for i in P.iterate(interval, step=1)])
    let_str = f"letting {name} be {{ {range_vals} }}\n"
    return let_str

def problem2essence(problem):
    essence = ""
    added_doms = []
    for lab, dom in problem.domains.items():
        dom_str = dom2essence(lab, dom)
        added_doms.append(lab)
        essence += dom_str
    sizes = problem.configuration.size.values
    if problem.configuration.size.values.upper == P.inf:
        ub = problem.universe.size() +1
        sizes = problem.configuration.size.values.replace(upper = ub)
    lenghts = P.iterate(sizes, step=1)
    sequence = problem.configuration.type == "sequence" 
    permutation =  problem.configuration.type == "permutation"
    subset = problem.configuration.type == "subset" 
    multiset =  problem.configuration.type == "multisubset"
    composition =  problem.configuration.type == "composition"
    partition =  problem.configuration.type == "partition"
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
                ub = "1" if permutation or subset else f"f_{uni}(e)"
                mset_constraint += f"\t\tsum([1 | i: int(1..l_{l}), {name}(i)=e]) <= {ub}\n "
                constraints.append(mset_constraint)
            for i in range(0,l):
                dom = ""
                for j, pf in enumerate(problem.pos_formulas):
                    if pf.pos-1 == i:
                        dlab = f"pf_{i}_{j}"
                        if dlab not in added_doms:
                            dom_str = dom2essence(dlab, pf.formula)
                            essence += dom_str
                            added_doms.append(dlab)
                        constraints.append(f"{name}({i}) in {dlab}")
            for i, cf in enumerate(problem.count_formulas):
                dlab = f"df_{i}"
                if dlab not in added_doms:
                    dom_str = dom2essence(dlab, cf.formula)
                    essence += dom_str
                    added_doms.append(dlab)
                range_vals = ",".join([str(i) for i in P.iterate(cf.values, step=1)])
                essence_l += f"letting vals_{i} be {{ {range_vals} }}\n"
                constraints.append(f"sum([1 | i: int(1..l_{l}), {name}(i) in {dlab}]) in vals_{i}")
            if sequence or permutation:
                essence_l += f"find {name} : sequence (size l_{l}) of {uni}\n"
            if multiset or subset:
                essence_l += f"find {name} : mset (size l_{l}) of {uni}\n"
            if len(constraints) >0:
                essence_l += "such that \n"
        else:
            myparts = [f"p{i}" for i in range(1,l+1)]
            essence_l += "letting myparts be new type enum {" + ",".join(myparts) + "}\n"
            essence_l += f"letting n be {problem.universe.size()}\n"
            nonempty = "forAll p: myparts.\n\t sum([put[e,p] | e:universe]) > 0"
            constraints.append(nonempty)
            alldistributed = "forAll e: universe.\n\t sum([put[e,p] | p: myparts]) = f_universe(e)"
            constraints.append(alldistributed)
            ub = problem.universe.size() - l + 1
            for i in range(1,l+1):
                dom = ""
                for j, pf in enumerate(problem.pos_formulas):
                    if pf.pos == i:
                        if pf.formula.size.values != P.closed(1,ub):
                            name = f"s_{i}"
                            essence_l += range2essence(pf.formula.size.values, name, ub)
                            size_constr = f"sum([put[e,p{i}] | e:{uni}]) in {name}"
                            constraints.append(size_constr)
                        else:
                            for k, cof in enumerate(pf.formula):
                                df = cof.formula
                                if dlab not in added_doms:
                                    dlab = f"df_{i}_{j}_{k}"
                                    dom_str = dom2essence(dlab, df)
                                    essence += dom_str
                                    added_doms.append(dlab)
                                range_name = f"vals_{i}_{j}_{k}"
                                essence_l += range2essence(cof.values, name)
                                constraints.append(f"sum([put[e,p{i}] | e:{dlab}]) in {range_name}")
            for i, cf in enumerate(problem.count_formulas):
                outer_range = range2essence(cf.values, f"vals_{i}_out", l)
                inner_range = range2essence(cf.formula.values, f"vals_{i}_in", ub)
                essence_l += outer_range + inner_range
                dlab = f"df_{i}"
                if dlab not in added_doms:
                    dom_str = dom2essence(dlab, cf.formula.formula)
                    essence += dom_str
                    added_doms.append(dlab)
                cf_constraint = f"|[p | p:myparts, sum([put[e,p] | e <- df_{i}]) in vals_{i}_in]| in vals_{i}_out"
                constraints.append(cf_constraint)

            essence_l += f"find put: matrix indexed by [universe, myparts] of int(0..n)\n"
            essence_l += "such that \n"
        constraint_str = "\n\t/\ ".join(constraints)
        essence_l += constraint_str
        essence_l += "\n"
        essence_l = essence + essence_l
        essence_lengths.append(essence_l)
    return essence_lengths


@timeout(TIMEOUT)
def run_asp(programs):
    n = 0
    for program in programs:
        # print(program)
        ctl = clingo.Control()
        ctl.configuration.solve.models = 0
        ctl.add("base", [],program)
        ctl.ground([("base", [])], context=Context())
        models = ctl.solve(yield_=True)
        n_length = sum(1 for _ in models)
        # print(n_length)
        n += n_length
    return n

@timeout(TIMEOUT)
def run_solver(problem):
    print("Running solver...")
    try:
        start = time.time()
        count = problem.solve(log=False)
        finish = time.time()
        print(f"Solver: {count} in {finish-start:.2f}s")
    except TimeoutError as e:
        print("CoSo timeout")

@timeout(TIMEOUT)
def run_sat(programs):
    n = 0
    gringo = os.path.join(asp_tools, "gringo")
    lp2normal = os.path.join(asp_tools, "lp2normal-2.18")
    lp2sat = os.path.join(asp_tools, "lp2sat-1.24")
    lp2atomic = os.path.join(asp_tools, "lp2atomic-1.17")
    input = os.path.join(asp_tools, "tmp.lp")
    out = os.path.join(asp_tools, "out.cnf")
    for program in programs:
        lp = open(input, "w+")
        lp.write(program)
        lp.close()
        p = Popen([f"{gringo} {input} | {lp2normal} | {lp2atomic} | {lp2sat} > {out}"], shell=True)
        p.wait()
        p = Popen([f"{sharp_sat} {out}"], shell=True, stdout=PIPE, stderr=PIPE)
        std_out, std_err = p.communicate()
        sol = std_out.decode('UTF-8')
        n_string = sol[sol.find("# solutions")+11:sol.find("# END")].replace('\n','').replace(' ','')
        n_program = int(n_string)
        # print(n_program)
        n += n_program
    os.remove(input)
    os.remove(out)
    return n

@timeout(TIMEOUT)
def run_essence_timeout(programs):
    n = 0
    for program in programs:
        # print(program)
        exec_conjure = os.path.join(conjure, "conjure")
        input = os.path.join(tools, "model.essence")
        model = open(input, "w+")
        model.write(program)
        model.close()
        conjure_env = os.environ.copy()
        conjure_env["PATH"] = os.path.abspath(conjure) + ":" + conjure_env["PATH"]
        p = Popen([f"{exec_conjure} solve -ac {input} --number-of-solutions=all --limit-time {TIMEOUT}"], shell=True, env=conjure_env)
        p.wait()
        solutions = [f for f in os.listdir(tools) if f.endswith(".solution")]
        n_prog = len(solutions)
        n += n_prog
        for s in solutions:
            os.remove(os.path.join(tools, s))
    return n

def run_essence(programs):
    n = run_essence_timeout(programs)
    for f in os.listdir(tools):
        if f.endswith(".solution"):
            os.remove(os.path.join(tools, f))
    return n
    
def compare(problem, name, sys_name, translate, run):
    print(f"Comparing on {name}")
    run_solver(problem)
    t = translate(problem)
    print(f"Running {sys_name}...")
    try:
        start = time.time()
        count = run(t)
        finish = time.time()
        print(f"{sys_name}: {count} in {finish-start:.2f}s")
    except TimeoutError as e:
        print(f"{sys_name} timeout")

# def problem2problog(problem):
#     problog = ":- use_module(library(aproblog)).\n\
#             :- use_semiring(\n\
#                 sr_plus,   % addition (arity 3)\n\
#                 sr_times,  % multiplication (arity 3)\n\
#                 sr_zero,   % neutral element of addition\n\
#                 sr_one,    % neutral element of multiplication\n\
#                 sr_pos,\n\
#                 sr_neg,    % negation of fact label\n\
#                 sr_negate, \n\
#                 true,      % requires solving disjoint sum problem?\n\
#                 true).    % requires solving neutral sum problem?\n\
#             sr_zero(0.0).\n\
#             sr_one(1.0).\n\
#             sr_plus(A, B, C) :- C is A + B.\n\
#             sr_times(A, B, C) :- C is A * B.\n\
#             sr_pos(A, B, A).\n\
#             sr_neg(A, B, 0).\n\
#             sr_negate(A, 1.0).\n"
#     n = problem.structure.size
#     u = problem.domains[problem.universe]
#     for i in range(1,n+1):
#         problog += ";".join(list(map(lambda e: f"1::v{i}({e})",portion.iterate(u.elements, step =1))))
#         problog += ".\n"
#     sequence = problem.structure.type == "sequence"
#     subset = problem.structure.type == "subset"
#     if sequence:
#         problog += f"sequences :- "
#     else:
#         problog += f"subsets :- "
#     var_pred = []
#     vars = []
#     for i in range(1,n+1):
#         var_pred.append(f"v{i}(V{i})")
#         vars.append(f"V{i}")
#     problog += ",".join(var_pred)
#     if sequence and problem.structure.spec or subset:
#         problog += ", "
#         if sequence and problem.structure.spec:
#             sign = "\="
#         else:
#             if problem.structure.spec:
#                 sign = "=<"
#             else:
#                 sign = "<"
#         ineq = []
#         for pair in itertools.product(vars,vars):
#             if pair[0] < pair[1]:
#                 ineq.append(f"{pair[0]} {sign} {pair[1]} ")
#         problog += ",".join(ineq)
#     problog += ".\n"
#     if sequence:
#         problog+="query(sequences).\n"
#     else:
#         problog+="query(subsets).\n"
#     return problog

# def compare2aproblog(problem ,name):
#     aproblog = problem2problog(problem)
#     probname = name[:-3] + "_prob.pl"
#     probpath = os.path.abspath(probname)
#     probfile = open(probname, "w")
#     probfile.write(aproblog)
#     probfile.close()
#     print("Running solver...")
#     start = time.time()
#     count = problem.solve(log=False)
#     finish = time.time()
#     print(f"Solver: {count} in {finish-start:.2f}s")
#     print("Running aProbLog...")
#     p = subprocess.Popen(
# 		["pyenv/bin/python3.8", "problog/problog-cli.py","-t 240", probname], 
# 		stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid)
#     finish = time.time()
#     try:
#         stdout, stderr = p.communicate()
#         finish = time.time()
#         out = stdout.decode('utf-8').split(":")
#         if out[0] == "Timeout exceeded\n" or out[0]=='':
#             print("\t TIMEOUT")
#         else:
#             out = " ".join(out[1].split())
#             probcount = int(float(out))
#             print(f"aProbLog: {probcount} in {finish-start:.2f}s")
#             if probcount == count:
#                 print("\t OK")
#             else:
#                 print("\t FAIL (check limitation of aProblog conversion)")
#     except subprocess.TimeoutExpired:
# 	    print("TIMEOUT")
# 	    p.terminate()
# 	    p.kill()
# 	    os.killpg(p.pid, signal.SIGINT)
# 	    print("KILLED")


def count_sols(filename):
    f = open(filename, "r")
    n = 0 
    for l in f.read():
        if l[-1] == ";":
            n +=1
    f.close()
    return n

# def compare2minizinc(problem, name):
#     minizinc = problem2minizinc(problem)
#     mininame = name[:-5] + "_mini.mzn"
#     minifile = open(mininame, "w")
#     minifile.write(minizinc)
#     minifile.close()
#     print("Running Solver...")
#     start = time.time()
#     sol = problem.solve(log=False)
#     count = sol.count
#     finish = time.time()
#     print(f"Solver: {count} in {finish-start:.2f}s")
#     print("Running Minizinc...")
#     start = time.time()
#     timeout = 240000
#     p = subprocess.Popen(
# 		["minizinc/bin/minizinc", "--time-limit", str(timeout), mininame, "-s", "-a"], 
# 		stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid)
#     try:
#         stdout, stderr = p.communicate()
#         finish = time.time()
#         if stderr != "":
#             lines = stdout.decode('utf-8').split("\n")
#             init_time = float(lines[-12].split('=')[1])
#             solve_time = float(lines[-11].split('=')[1])
#             tot_time = init_time + solve_time
#             minicount = int(lines[-10].split('=')[1])
#         else:
#             minicount = 0
#             tot_time = start - finish
#         print(f"Minizinc: {minicount} in {tot_time:.2f}s")
#         if minicount == count and finish-start < (timeout-0.1*timeout)/1000:
#             print("\t OK")
#         else:
#             if finish-start > (timeout-0.1*timeout)/1000:
#                 print("\t TIMEOUT")
#             else:
#                 print("\t FAIL")
#     except subprocess.TimeoutExpired:
# 	    print("FAIL")
#         p.terminate()
#         p.kill()
#         os.killpg(p.pid, signal.SIGINT)
#         print("KILLED")

# def generate_unconstrained():
#     size_domain = [4,9,14]
#     size_struct = [4,7,10]
#     cases = [("sequence", True, "false"), ("permutation", True, "true"), ("subset", False, "false")]
#     for d in size_domain:
#         for s in size_struct:
#             if d>s:
#                 for c in cases:
#                     name, seq, spec = c 
#                     print(f"{name}, size domain={d}, size structure={s}")
#                     problem = generate_problem(3,d,seq,s,spec,False,False)
#                     print(problem)
#                     Path("tests/unconstrained").mkdir(parents=True, exist_ok=True)
#                     filename = f"tests/unconstrained/{name}_{d}_{s}.test"
#                     file = open(filename, "w")
#                     file.write(problem)
#                     file.close()

def generate_constrained(folder, pconstr, cconstr):
    size_domain = [10,15,20]
    size_struct = [5,10,15]
    cases = ["sequence", "permutation", "subset", "multisubset", "partition", "composition"]
    for d in size_domain:
        domains_upto = d // 3
        for s in size_struct:
            for case in cases:
                case_path = os.path.join(folder,case)
                os.makedirs(case_path, exist_ok=True)
                if d>s:
                    # print(f"{name}, size domain={d}, size structure={s}")
                    problem = generate_problem(case, domains_upto, d, s, pconstr, cconstr)
                    filename = f"{case}_{d}_{s}.test"
                    path = os.path.join(case_path, filename)
                    file = open(path, "w")
                    file.write(problem)
                    file.close()
                    parser = Parser(str(path))
                    parser.parse()
                    problem_asp = problem2asp(parser.problem)
                    filename_asp = f"{case}_{d}_{s}.lp"
                    path_asp = os.path.join(case_path, filename_asp)
                    file_asp = open(path_asp, "w")
                    file_asp.write(problem_asp)
                    file_asp.close()

# def test_folder(folder, problog, from_file = ""):
#     for filename in os.listdir(folder):
#         if filename.endswith(".test"):
#             if from_file =="" or from_file !="" and filename>from_file:
#                 print(f"Test {filename}:")
#                 parser = Parser(os.path.join(folder,filename))
#                 problem = parser.parsed
#                 compare2minizinc(problem)
#                 if problog:
#                     compare2aproblog(problem)

def test_folder(folder, asp, sat, essence):
    for filename in os.listdir(folder):
        if filename.endswith(".test") or filename.endswith(".pl"):
            print(f"Test {filename}:")
            parser = Parser(os.path.join(folder,filename))
            try:
                parser.parse()
                problem = parser.problem
                name = os.path.join(folder,filename)
                if essence:
                    compare(parser.problem, name, "Essence", problem2essence, run_essence)
                if asp:
                    compare(parser.problem, name, "Clingo", problem2asp, run_asp)
                if sat:
                    compare(parser.problem, name, "SharpSAT", problem2asp, run_sat)
                if not essence and not asp and not sat:
                    run_solver(problem)
            except EmptyException:
                print("Empty: skipping")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Run solver on file')
    parser.add_argument('-l', action='store_true', help='Logging')
    parser.add_argument('-g', help='Geneate random tests in given folder')
    parser.add_argument('-m', default="minizinc", help='Minizinc directory')  
    parser.add_argument('-t', "--timeout", type=int, default=TIMEOUT, help="Set timeout in seconds")  
    parser.add_argument('--test-folder', help='Run tool comparison on files in folder')
    parser.add_argument('--minizinc', action='store_true', help="Compare with 'minizinc'")
    parser.add_argument('--asp', action='store_true', help="Compare with 'asp'")
    parser.add_argument('--sat', action='store_true', help="Compare with 'sharpSAT'")
    parser.add_argument('--essence', action='store_true', help="Compare with 'essence'")
    parser.add_argument('--noposconstr', action='store_false', help="Disable positional constraint generation when -g")
    parser.add_argument('--nocountconstr', action='store_false', help="Disable counting constraint generation when -g")
    args = parser.parse_args()
    # print(args)
    if args.f:
        parser = Parser(args.f)
        parser.parse()
        # print(parser.problem)
        # if args.minizinc:
        #     compare2minizinc(parser.problem, args.f)
        if args.asp:
            compare(parser.problem, args.f, "Clingo", problem2asp, run_asp)
        if args.sat:
            compare(parser.problem, args.f, "SharpSAT", problem2asp, run_sat)
        if args.essence:
            compare(parser.problem, args.f, "Essence", problem2essence, run_essence)
        if not args.asp and not args.sat and not args.essence:
            sol = parser.problem.solve(log=args.l)
            print(f"Count: {sol}")
    elif args.g:
        generate_constrained(args.g, args.noposconstr, args.nocountconstr)
    elif args.test_folder:
        # ap = 'aproblog' in args.compare
        test_folder(args.test_folder, args.asp, args.sat, args.essence)
    else:
        pass
