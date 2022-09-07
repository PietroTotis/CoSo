import portion as Int
from solver import Solver
from level_1 import SetFormula
from configuration import CSize
from sharpCSP import Solution
from util import *


class EmptyException(Exception):
    def __init__(self, message):
        super().__init__(message)


class Problem(object):
    """
    Represents a problem as a collection of domains, a target configuration, a set of choice formulas and a set of constraint formulas.

    Attributes
    ----------
    universe : str
        the name of the largest domain that includes all the others (assumes it was declared)
    entity_map : dict str->int
        maps each entity/constant to an integer (manipulated as intervals)
    """

    def __init__(self):
        self.universe = None
        # self.agg_formulas = []
        self.pos_formulas = []
        self.count_formulas = []
        self.domains = {}
        self.entity_map = {}
        self.configuration = None
        self.internal_copies = {}
        self.label_map = {}

    # def add_choice_formula(self, head, body):
    #     struct_name = head.args[0]
    #     formula = head.args[-1]
    #     df = self.compute_dom(formula)
    #     if head.functor == "pos":
    #         pos = head.args[1]
    #         c = PosFormula(struct_name, pos, df)
    #     if head.functor == "in":
    #         c = InFormula(self.configuration.name, df)
    #     self.choice_formulas.append(c)

    def add_domain(self, dom):
        self.domains[dom.name] = dom

    def add_pos_formula(self, chf):
        self.pos_formulas.append(chf)

    def add_agg_formula(self, chf):
        self.agg_formulas.append(chf)

    def add_counting_formula(self, cof):
        # cformula = self.build_cof(cof)
        subsets = self.configuration.type in ["partition", "composition"]
        count_prop = isinstance(cof.formula, SetFormula)
        # ignore recursive parsing of counting formulas
        if not (subsets and count_prop):
            self.count_formulas.append(cof)

    def add_entity(self, e):
        """
        Map an entity e to its internal representation (id)

        Args:
            e (str): the label of an entity

        Returns:
            int: the id representing the entity internally
        """
        name = str(e)
        if name in self.entity_map:
            return self.entity_map[name]
        else:
            i = len(self.entity_map) + 1
            self.entity_map[name] = i
            copy_ids = name.split("__")
            if len(copy_ids) > 1:
                base_id, _ = copy_ids
                if base_id in self.internal_copies:
                    self.internal_copies[base_id] += [i]
                else:
                    self.internal_copies[base_id] = [i]
            return i

    def get_entity(self, e):
        name = str(e)
        base_id = self.entity_map.get(name)
        if base_id is None:
            raise Exception(f"Unknown label {name}")
        copy_ids = self.internal_copies.get(name, [])
        return [base_id] + copy_ids

    # def add_query(self, q):
    #     self.queries.append(q)

    # def add_size(self, sf):
    #     name = sf.args[0]
    #     op = sf.args[1]
    #     formula = sf.args[2]
    #     s = self.build_size(name,op,formula)
    #     self.configuration.size = s

    # def add_configuration(self, struct):
    #     self.configuration = struct

    # def build_cof(self, cof):
    #     struct_name = cof.args[0]
    #     problog_formula = cof.args[1]
    #     op = cof.args[2].functor
    #     val = cof.args[3].compute_value()
    #     formula = self.compute_formula(problog_formula)
    #     interval = self.get_interval(op, val)
    #     return CountingFormula(formula, interval)

    # def build_size(self, name, op, val):
    #     n = val.compute_value()
    #     op = op.functor
    #     interval = self.get_interval(op,n)
    #     return SizeFormula(name, interval)

    def compute_dom(self, sformula):
        """
        Given a set-formula expands and computes the corresponding domain
        Parameters
        ----------
        sformula : a set formula : an And/Or/Not (possibly nested)
                  or a string corresponding to one of the domains
        """
        if isinstance(sformula, And):
            lf = self.compute_dom(sformula.left)
            rf = self.compute_dom(sformula.right)
            domain = lf & rf
        elif isinstance(sformula, Or):
            lf = self.compute_dom(sformula.left)
            rf = self.compute_dom(sformula.right)
            domain = lf | rf
        elif isinstance(sformula, Not):
            arg = self.compute_dom(sformula.child)
            domain = arg.neg()
        else:  # base cases: universe, arrangement, user-defined, single element
            if sformula in self.domains:
                domain = self.domains[sformula]
            elif (
                sformula == "universe"
                or sformula == "part"
                or self.configuration is not None
                and sformula == self.configuration.name
            ):
                domain = self.universe
            else:
                id = self.get_entity(sformula)
                if id is None:
                    id = ""  # sformula # dummy
                    single = Int.closed(0, P.inf)
                    # raise Exception(f"Unknown constant {sformula}")
                else:
                    single = Int.singleton(id)
                dist = Int.IntervalDict()
                dist[single] = True
                domain = SetFormula(id, dist, self.universe)
        return domain

    def compute_universe(self):
        """
        Once all sets are added, compute the union and update the domain formulas
        with all info about the universe
        """
        dom_iter = iter(self.domains.values())
        universe = next(dom_iter)
        for d in dom_iter:
            universe = universe | d
        universe.name = "universe"
        dom_iter = iter(self.domains.values())
        self.label_map = {v: k for k, v in self.entity_map.items()}
        universe.set_labels(self.label_map)
        for d in dom_iter:
            d.set_universe(universe)
            d.set_labels(self.label_map)
        for pf in self.pos_formulas:
            pf.set_universe(universe)
            pf.set_labels(self.label_map)
        for cof in self.count_formulas:
            cof.set_universe(universe)
            cof.set_labels(self.label_map)
            # cof.formula.universe = universe
        self.universe = universe
        self.universe.set_universe(universe)
        print("dsdsmpc", self.universe, self.universe.universe)

    # def compute_formula(self, formula):
    #     if formula.functor == "size":
    #         name = formula.args[0]
    #         op = formula.args[1]
    #         n = formula.args[2]
    #         return self.build_size(name,op,n)
    #     elif formula.functor == "count":
    #         return self.build_cof(formula)
    #     else:
    #         return self.compute_dom(formula)

    # def get_entity(self, e):
    #     if e in self.entity_map:
    #         return self.entity_map[e]
    #     else:
    #         return None

    def get_interval(self, op, n):
        """
        Convert a comparison operator and a number to interval

        Args:
            op (str): a string describing a comparison
            n (int): reference value

        Returns:
            Int: interval of values
        """
        if op == "<":
            interval = Int.closedopen(0, n)
        elif op == "<=":
            interval = Int.closed(0, n)
        elif op == ">":
            interval = Int.open(n, Int.inf)
        elif op == ">=":
            interval = Int.closedopen(n, Int.inf)
        elif op == "=":
            interval = Int.singleton(n)
        else:
            interval = Int.closedopen(0, n) | Int.open(n, Int.inf)
        return interval

    def solve(self, log=True):
        """
        Do some sanity checks on the configuration and then call the solver

        Args:
            log (bool, optional): enable logging. Defaults to True.

        Returns:
            Solution: solution
        """
        if self.configuration is None or len(self.domains) == 0:
            raise EmptyException("Empty problem!")
        else:
            if self.configuration.size is None:
                vals = Int.closed(1, self.universe.size())
                self.configuration.size = CSize("unconstrained", vals)
            s = Solver(self)
            return s.solve()

    def __str__(self):
        s = ""
        for d in self.domains.values():
            s += f"{d.name}: {d}\n"
        for cf in self.pos_formulas:
            s += f"{cf}\n"
        for f in self.count_formulas:
            s += f"{f}\n"
        for f in self.agg_formulas:
            s += f"{f}\n"
        s += f"{self.configuration}\n"
        return s
