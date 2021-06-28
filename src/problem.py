import portion as P
from configuration import *
from solver import Solver
from formulas import *
from util import *
 

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
        self.agg_formulas = []
        self.pos_formulas = []
        self.count_formulas = []
        self.domains = {}
        self.entity_map = {}
        self.configuration = None

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
        subsets = self.configuration.type in ["partition","composition"]
        count_prop = isinstance(cof.formula, DomainFormula)
        # ignore recursive parsing of counting formulas
        if not (subsets and count_prop):
            self.count_formulas.append(cof)

    def add_entity(self, e):
        if str(e) in self.entity_map:
            return self.entity_map[str(e)]
        else:
            i = len(self.entity_map) +1
            self.entity_map[str(e)] = i
            return i

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
            domain =  arg.neg()
        else: # base cases: universe, arrangement, user-defined, single element
            if sformula in self.domains:
                domain = self.domains[sformula]
            elif sformula == "universe" or self.configuration is not None and sformula == self.configuration.name:
                domain = self.universe
            else:
                id = self.get_entity(sformula)
                if id is None:
                    id = sformula # dummy
                    single = portion.closed(0,P.inf)
                    # raise Exception(f"Unknown constant {sformula}")
                else:
                    single = portion.singleton(id)
                dist = portion.IntervalDict()
                dist[single] = True
                domain = DomainFormula(id, dist, self.universe)
        return domain

    def compute_universe(self):
        dom_iter = iter(self.domains.values())
        universe = next(dom_iter)
        for d in dom_iter:
            universe = universe | d
        # universe.name = "universe"
        dom_iter = iter(self.domains.values())
        for d in dom_iter:
            d.universe = universe
        # if self.configuration.df is None:
        #     self.configuration.df  = universe
        for pf in self.pos_formulas:
            pf.formula.universe = universe
        for cof in self.count_formulas:
            cof.formula.universe = universe
        self.universe = universe
        self.universe.universe = universe

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

    def get_entity(self, e):
        if e in self.entity_map:
            return self.entity_map[e]
        else:
            return None

    def get_interval(self,op,n):
        if op == "<":
            interval = portion.closedopen(0,n)
        elif op == "<=":
            interval = portion.closed(0,n)
        elif op == ">":
            interval = portion.open(n,portion.inf)
        elif op == ">=":
            interval = portion.closedopen(n,portion.inf)
        elif op == "=":
            interval = portion.singleton(n)
        else:
            interval = portion.closedopen(1,n) | portion.open(n,portion.inf)
        return interval

    def solve(self, log=True):
        if self.configuration is None or len(self.domains) == 0:
            print("empty problem!")
            return 0
        else:
            if self.configuration.size is None:
                vals = portion.closed(1, self.universe.size())
                self.configuration.size = SizeFormula("unconstrained", vals)
            s = Solver(self)
            return s.solve(log)

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
