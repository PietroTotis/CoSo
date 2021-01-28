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
        self.universe = problem.structure.df.domain
        self.size  = problem.structure.size
        if self.size.values.upper == P.inf:
            n = self.universe.size()
            self.size.values = self.size.values.replace(upper=n, right=P.CLOSED)
        self.type = problem.structure.type

    def solve(self, log=True):
        count = Solution(0,[])
        var_dom = DomainFormula(self.universe, "universe", self.universe)
        for n in self.size:
            if self.type in ["sequence", "subset"]:
                vars = [var_dom]*n
                csp = SharpCSP(vars, self.type, self.problem.choice_formulas, self.problem.count_formulas, self.problem.structure.spec, var_dom)
                count += csp.solve(log)
            else:
                count = 0
                ub_size = self.universe.size() - n + 1
                size = SizeFormula("universe", P.closed(1,ub_size))
                vars = [LiftedSet(f"part. of {var_dom}", size)]*n
            csp = SharpCSP(vars, self.type, self.problem.choice_formulas, self.problem.count_formulas, self.problem.structure.spec, var_dom)
            count += csp.solve(log)
        return count
