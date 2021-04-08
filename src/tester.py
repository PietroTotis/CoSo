import os
import argparse
import time
import portion
import subprocess
import pyinotify
import signal
import random
import itertools
from pathlib import Path
from parser import Parser
from formulas import PosFormula, And, Or, Not

ops = [">","<","<=",">=","!=","="]

class ModHandler(pyinotify.ProcessEvent):
    def process_IN_CLOSE_WRITE(self, evt):
        pass

def get_random_complex_dom(n_domains):
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

def generate_problem(n_domains,u_size,is_sequence,struct_size,spec,choice_constraints,counting_constraints):
    random.seed(4) #just less random problems with 0 solutions
    problem = ""
    universe_size = u_size
    problem += f"uni = [1:{universe_size}];\n"
    for i in range(1,n_domains):
        base = random.randint(1,universe_size-1)
        upper = random.randint(base+1,universe_size)
        problem+= f"dom{i} = [{base}:{upper}];\n"
    if is_sequence:
        problem += f"a in [{spec} uni];\n"
    else:
        problem += f"a in {{{spec} uni}};\n"
    problem += f"#a = {struct_size};\n"
    if choice_constraints and is_sequence:
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
            problem += f"#{dom} {op} {n};\n"
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

def get_n_vars(n):
    vars = []
    for c in range(ord('A'),ord('Z')):
        if len(vars) < n:
            vars.append(c)


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

def compare2minizinc(problem, name):
    minizinc = problem2minizinc(problem)
    mininame = name[:-5] + "_mini.mzn"
    minifile = open(mininame, "w")
    minifile.write(minizinc)
    minifile.close()
    print("Running Solver...")
    start = time.time()
    sol = problem.solve(log=False)
    count = sol.count
    finish = time.time()
    print(f"Solver: {count} in {finish-start:.2f}s")
    print("Running Minizinc...")
    start = time.time()
    timeout = 240000
    p = subprocess.Popen(
		["minizinc/bin/minizinc", "--time-limit", str(timeout), mininame, "-s", "-a"], 
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid)
    try:
        stdout, stderr = p.communicate()
        finish = time.time()
        if stderr != "":
            lines = stdout.decode('utf-8').split("\n")
            init_time = float(lines[-12].split('=')[1])
            solve_time = float(lines[-11].split('=')[1])
            tot_time = init_time + solve_time
            minicount = int(lines[-10].split('=')[1])
        else:
            minicount = 0
            tot_time = start - finish
        print(f"Minizinc: {minicount} in {tot_time:.2f}s")
        if minicount == count and finish-start < (timeout-0.1*timeout)/1000:
            print("\t OK")
        else:
            if finish-start > (timeout-0.1*timeout)/1000:
                print("\t TIMEOUT")
            else:
                print("\t FAIL")
    except subprocess.TimeoutExpired:
	    print("FAIL")
	    p.terminate()
	    p.kill()
	    os.killpg(p.pid, signal.SIGINT)
	    print("KILLED")

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

def generate_constrained(folder):
    size_domain = [10,15,20]
    size_struct = [5,10,15]
    cases = [("sequence", True, "||"), ("permutation", True, "|"), ("subset", False, "|")]
    for d in size_domain:
        for s in size_struct:
            for c in cases:
                if d>s:
                    name, seq, spec = c 
                    print(f"{name}, size domain={d}, size structure={s}")
                    problem = generate_problem(3,d,seq,s,spec,True,True)
                    print(problem)
                    filename = f"{name}_{d}_{s}.test"
                    path = os.path.join(folder,filename)
                    file = open(path, "w")
                    file.write(problem)
                    file.close()

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

def test_folder(folder, aproblog, minizinc):
    for filename in os.listdir(folder):
        if filename.endswith(".test") or filename.endswith(".pl"):
            print(f"Test {filename}:")
            parser = Parser(os.path.join(folder,filename))
            parser.parse()
            problem = parser.problem
            if minizinc:
                compare2minizinc(problem, os.path.join(folder,filename))
            if aproblog:
                compare2aproblog(problem,  os.path.join(folder,filename))
            if not minizinc and not aproblog:
                sol = problem.solve(log=False)
                print(f"Count: {sol}") 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Run solver on file')
    parser.add_argument('-g', help='Geneate random tests in given folder')
    parser.add_argument('-m', default="minizinc", help='Minizinc directory')    
    parser.add_argument('--test-folder', help='Run tool comparison on files in folder')
    parser.add_argument('--minizinc', action='store_true', help="Compare with 'minizinc'")
    args = parser.parse_args()
    if args.f:
        parser = Parser(args.f)
        parser.parse()
        print(parser.problem)
        if args.minizinc:
            compare2minizinc(parser.problem, args.f)
        else:
            sol = parser.problem.solve(log=False)
            print(f"Count: {sol}")
    elif args.g:
        generate_constrained(args.g)
    elif args.test_folder:
        # ap = 'aproblog' in args.compare
        test_folder(args.test_folder, False,  args.minizinc)
    else:
        pass
