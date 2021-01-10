import portion as P
import math 
import itertools

from problog.logic import Constant

from formulas import *
from sharpCSP import SharpCSP, Solution
from structure import Domain, LiftedSet

class Solver(object):
    """
    Sets up the decision variables from the problem
    """
    def __init__(self,problem):
        self.problem = problem
        self.universe = problem.universe
        self.size = problem.structure.size
        self.type = problem.structure.type

    def solve(self, log=True):
        count = Solution(0,[])
        var_dom = DomainFormula(self.universe.name, Term("universe"), self.universe)
        for n in self.size:
            if self.type in ["sequence", "subset"]:
                vars = [var_dom]*n
            else:
                ub_size = self.universe.size() - n + 1
                size = SizeFormula("universe", P.closed(1,ub_size))
                vars = [LiftedSet(f"part. of {var_dom}", size)]*n
            csp = SharpCSP(vars, self.type, self.problem.choice_formulas, self.problem.count_formulas, self.problem.structure.spec, var_dom)
            count += csp.solve(log)
        return count
