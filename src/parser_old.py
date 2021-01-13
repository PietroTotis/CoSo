import argparse
import portion
import operator, functools

from problog.parser import PrologParser
from problog.parser import ParseError
from problog.program import PrologFile
from problog.logic import Clause, term2list

from problem import *
from formulas import *

reserved = ["count", "query", "size", "pos", "in", "part"]

class Parser(object):

    def __init__(self, file):
        self.parsed = Problem()
        self.parse_file(file)

    def add_clause(self, cl):
        self.parsed.add_choice_formula(cl.head,cl.body)

    def add_domain(self, stmt):
        """ 
        Adds to the problem a domain expressed ether as an interval 'name((in)distinguishable,[lb,ub])' 
        or as an enumeration 'name((in)distinguishable,[e1,...en])'
        
        Parameters
        ----------
        problem : Problem to which the statement belongs
        stmt : a ProbLog clause
        """
        if len(stmt.args) == 1:
            distinguishable = True
            elems = term2list(stmt.args[0])
        else:
            distinguishable = stmt.args[0] == "distinguishable"
            elems = term2list(stmt.args[1])
        if len(elems) == 2 and isinstance(elems[0],int) and isinstance(elems[1],int):
            ivs = [portion.closed(elems[0], elems[1])]
        else:
            entities = list(map(self.parsed.add_entity, elems))
            entities.sort()
            ivs = []
            i = 0
            while i< len(entities):
                low = entities[i]
                while i <len(entities)-1 and entities[i]+1 == entities[i+1]: i +=1
                hi = entities[i]
                if hi - low >=1:
                    ivs.append(portion.closed(low, hi))
                else:
                    ivs.append(portion.singleton(low))
                i += 1
        elements = functools.reduce(lambda a,b: a.union(b), ivs, portion.empty())
        dist = P.IntervalDict()
        dist[elements] = distinguishable
        d = Domain(stmt.functor, elements, dist)
        self.parsed.add_domain(d)


    def add_statement(self, stmt):
        if stmt.functor == 'structure':
            s = Structure(*stmt.args)
            self.parsed.add_structure(s)
        elif stmt.functor not in reserved:
            self.add_domain(stmt)
        elif stmt.functor == "count":
            self.parsed.add_counting_formula(stmt)
        elif stmt.functor == "query":
            q = stmt.args[0]
            self.parsed.add_query(q)
        elif stmt.functor == "size":
            self.parsed.add_size(stmt)
        elif stmt.functor in ["pos", "in", "part"]:
            self.parsed.add_choice_formula(stmt, None)
        else:
            pass

    def parse_file(self, filename):
        program = PrologFile(filename)
        for stmt in program:
            if isinstance(stmt, Clause):
                self.add_clause(stmt)
            else:
                self.add_statement(stmt)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('program', help='file name')
    args = parser.parse_args()
    p = Parser(args.program)
    print(p.parsed)
    sol = p.parsed.solve()
    print("count: ", sol)