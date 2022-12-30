import portion as P

from configuration import *
from count import Zero
from sharpCSP import SharpCSP
from level_2 import LiftedSet
from logger import ProblemLog


class Solver(object):
    """
    Sets up the decision variables from the problem
    """

    def __init__(self, problem, debug=False):
        self.debug = debug
        self.problem = problem
        self.universe = problem.universe
        self.size = problem.configuration.size
        if self.size.values.upper == P.inf:
            n = self.universe.size()
            self.size.values = self.size.values.replace(upper=n, right=P.CLOSED)
        self.type = problem.configuration.type
        self.log = ProblemLog(
            universe=self.universe,
            debug=debug,
            type=self.type,
            pos_constraints=self.problem.pos_constraints,
            constraints=problem.constraints,
            configuration=problem.configuration,
        )
        self.log.description("Root problem")

    def solve(self):
        if not self.trivial_unsat():
            count = Zero(self.log)
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
                    self.problem.pos_constraints,
                    self.problem.constraints,
                    self.universe,
                    caption=f"Configuration of size {n}",
                    lvl=1,
                    debug=self.debug,
                )
                size_count = csp.solve()
                self.log.add_subproblem("add", size_count.log)
                count += size_count
            self.log.solution = count
        else:
            count = Zero(self.log)
        return count

    def trivial_unsat(self):
        """
        Check if there is some constraint unsatisfiable because of domains, regardless of the number of variables

        Args:
            log (bool, optional): Logging. Defaults to True.
        """
        if self.type in ["sequence", "subset", "permutation", "multisubset"]:
            for count_constr in self.problem.constraints:
                lower_b = count_constr.values.lower
                atleast = (
                    lower_b if count_constr.values.left == P.CLOSED else lower_b + 1
                )
                available_elems = count_constr.formula.size()
                if available_elems < atleast:
                    self.log.description = f"Trivially unsat because we want at least {atleast} {count_constr.formula} but we have only {available_elems}"
                    return True
        else:
            # TODO
            return False
