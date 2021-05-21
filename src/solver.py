import portion as P
import math 
import itertools

from formulas import *
from sharpCSP import SharpCSP, Solution
from configuration import Domain, LiftedSet

class Solver(object):
    """
    Sets up the decision variables from the problem
    """
    def __init__(self,problem):
        self.problem = problem
        self.universe = problem.configuration.df
        self.size  = problem.configuration.size
        if self.size.values.upper == P.inf:
            n = self.universe.size()
            self.size.values = self.size.values.replace(upper=n, right=P.CLOSED)
        self.type = problem.configuration.type

    def solve(self, log=True):
        count = Solution(0,[])
        for n in self.size:
            if self.type in ["sequence", "subset", "permutation", "multisubset"]:
                vars = [self.universe]*n
            else:
                ub_size = self.universe.size() - n + 1
                size = SizeFormula("universe", P.closed(1,ub_size))
                vars = [LiftedSet(self.universe, size)]*n
            csp = SharpCSP(vars, self.type, self.problem.pos_formulas, self.problem.count_formulas, self.universe)
            count += csp.solve(log)
        return count
