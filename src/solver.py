import portion as P

from .configuration import *
from .count import Zero, Solution
from .sharpCSP import SharpCSP
from .level_2 import LiftedSet
from .logger import ProblemLog


class Solver(object):
    """
    Sets up the decision variables from the problem
    """

    def __init__(self, problem, debug=False):
        self.debug = debug
        self.problem = problem
        self.universe = problem.universe
        n = self.universe.size()
        problem.configuration.size.update_upper_bound(n)
        for c in self.problem.constraints:
            c.update_upper_bound(n)
        self.size = problem.configuration.size
        self.config = problem.configuration
        self.log = ProblemLog(
            universe=self.universe,
            debug=debug,
            config=self.config,
            pos_constraints=self.problem.pos_constraints,
            constraints=problem.constraints,
            configuration=problem.configuration,
        )
        self.log.description("Root problem")

    def solve(self):
        unsat, msg = self.trivial_unsat()
        if not unsat:
            count = Zero()
            for n in self.size:
                if self.config.lvl1():
                    vars = [self.universe] * n
                else:
                    ub_size = self.universe.size() - n + 1
                    size = CSize("universe", P.closed(1, ub_size))
                    vars = [LiftedSet(self.universe, size)] * n
                csp = SharpCSP(
                    vars,
                    self.config,
                    self.problem.pos_constraints,
                    self.problem.constraints,
                    self.universe,
                    caption=f"Configuration of size {n}",
                    lvl=1,
                    debug=self.debug,
                )
                size_sol = csp.solve()

                count += size_sol.count
                self.log.add_subproblem("add", size_sol.log)
        else:
            count = Zero(tip=msg)

        return Solution(count, self.log)

    def trivial_unsat(self):
        """
        Check if there is some constraint unsatisfiable because of domains, regardless of the number of variables

        Args:
            log (bool, optional): Logging. Defaults to True.
        """
        if self.config.lvl1():
            for count_constr in self.problem.constraints:
                lower_b = count_constr.values.lower
                atleast = (
                    lower_b if count_constr.values.left == P.CLOSED else lower_b + 1
                )
                available_elems = count_constr.formula.size()
                if available_elems < atleast:
                    msg = f"Trivially unsat because we want at least {atleast} {count_constr.formula} but we have only {available_elems}"
                    self.log.description = msg
                    return (True, msg)
            return (False, "")
        else:
            # TODO
            return (False, "")
