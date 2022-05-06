import portion as P

from configuration import *
from sharpCSP import SharpCSP, Solution
from level_2 import LiftedSet


class Solver(object):
    """
    Sets up the decision variables from the problem
    """

    def __init__(self, problem):
        self.problem = problem
        self.universe = problem.configuration.df
        self.size = problem.configuration.size
        if self.size.values.upper == P.inf:
            n = self.universe.size()
            self.size.values = self.size.values.replace(upper=n, right=P.CLOSED)
        self.type = problem.configuration.type

    def solve(self, log=True):
        if not self.trivial_unsat():
            count = Solution(0, [], subproblems=0)
            for n in self.size:
                if self.type in ["sequence", "subset", "permutation", "multisubset"]:
                    vars = [self.universe] * n
                else:
                    ub_size = self.universe.size() - n + 1
                    size = CSize("universe", P.closed(1, ub_size))
                    vars = [LiftedSet(self.universe, size)] * n
                csp = SharpCSP(
                    vars,
                    self.type,
                    self.problem.pos_formulas,
                    self.problem.count_formulas,
                    self.universe,
                )
                count += csp.solve(log)
        else:
            count = Solution(0, [], subproblems=1)
        return count

    def trivial_unsat(self, log=True):
        """
        Check if there is some constraint unsatisfiable because of domains, regardless of the number of variables

        Args:
            log (bool, optional): Logging. Defaults to True.
        """
        if self.type in ["sequence", "subset", "permutation", "multisubset"]:
            for count_constr in self.problem.count_formulas:
                lower_b = count_constr.values.lower
                atleast = (
                    lower_b if count_constr.values.left == P.CLOSED else lower_b + 1
                )
                available_elems = count_constr.formula.size()
                if available_elems < atleast:
                    if log:
                        print(
                            f"\tTrivially unsat because we want at least {atleast} {count_constr.formula} but we have only {available_elems}"
                        )
                    return True
        else:
            # TODO
            return False
