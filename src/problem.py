
from structure import *
from solver import Solver
from formulas import *
 

class Problem(object):
    """
    Represents a problem as a collection of domains, a target structure, a set of choice formulas and a set of constraint formulas (queries?).

    Attributes
    ----------
    universe : str
        the name of the largest domain that includes all the others (assumes it was declared)
    entity_map : dict str->int
        maps each entity/constant to an integer (manipulated as intervals)
    """
    
    def __init__(self):
        self.choice_formulas = []
        self.universe = None
        self.count_formulas = []
        self.domains = {}
        self.entity_map = {}
        self.queries = []
        self.structure = None

    def add_choice_formula(self, head, body):
        struct_name = head.args[0]
        formula = head.args[-1]
        df = self.compute_dom(formula)
        if head.functor == "pos":
            pos = head.args[1]
            c = PosFormula(struct_name, pos, df)
        if head.functor == "in":
            c = InFormula(self.structure.name, df)
        self.choice_formulas.append(c)

    def add_domain(self, dom):
        self.domains[dom.name] = dom
        if self.universe is None:
            self.universe = dom.name
        else: 
            universe = self.domains[self.universe]
            if universe in dom:
                self.universe = dom.name

    def add_counting_formula(self, cof):
        cformula = self.build_cof(cof)
        self.count_formulas.append(cformula)

    def add_entity(self, e):
        if str(e) in self.entity_map:
            return self.entity_map[str(e)]
        else:
            i = len(self.entity_map) +1
            self.entity_map[str(e)] = i
            return i

    def add_query(self, q):
        self.queries.append(q)

    def add_size(self, sf):
        name = sf.args[0]
        op = sf.args[1]
        formula = sf.args[2]
        s = self.build_size(name,op,formula)
        self.structure.size = s

    def add_structure(self, struct):
        self.structure = struct

    def build_cof(self, cof):
        struct_name = cof.args[0]
        problog_formula = cof.args[1]
        op = cof.args[2].functor
        val = cof.args[3].compute_value()
        formula = self.compute_formula(problog_formula)
        interval = self.get_interval(op, val)
        return CountingFormula(formula, interval)
    
    def build_size(self, name, op, val):
        n = val.compute_value()
        op = op.functor
        interval = self.get_interval(op,n)
        return SizeFormula(name, interval)

    def compute_dom(self, dformula):
        """
        Given a formula expands and computes the corresponding domain
        Parameters
        ----------
        dformula : a domain formula : a ProbLog predicate inter/union/not (possibly nested)
                  or a ProbLog Term with functor one of the domains
        """
        cont = self.domains[self.universe]
        if dformula.functor == "inter":
            lf = self.compute_dom(dformula.args[0])
            rf = self.compute_dom(dformula.args[1])
            domain = lf.domain & rf.domain
        elif dformula.functor == "union":
            lf = self.compute_dom(dformula.args[0])
            rf = self.compute_dom(dformula.args[1])
            domain = lf.domain | rf.domain
        elif dformula.functor == "not":
            arg = self.compute_dom(dformula.args[0])
            domain =  arg.neg().domain
        else:
            if dformula.functor in self.domains:
                domain = self.domains[str(dformula)]
            elif dformula.functor == "partition":
                domain = Domain("", portion.empty())
            else:
                id = self.get_entity(dformula.functor)
                if id is None:
                    raise Exception(f"Unknown constant {dformula.functor}")
                domain = Domain(id, portion.singleton(id))
        df = DomainFormula(cont, dformula, domain)
        return df

    def compute_formula(self, formula):
        if formula.functor == "size":
            name = formula.args[0]
            op = formula.args[1]
            n = formula.args[2]
            return self.build_size(name,op,n)
        elif formula.functor == "count":
            return self.build_cof(formula)
        else:
            return self.compute_dom(formula)

    def get_entity(self, e):
        if e in self.entity_map:
            return self.entity_map[e]
        else:
            return None

    def get_interval(self,op,n):
        if op == "<":
            interval = portion.closedopen(0,n)
        elif op == "=<":
            interval = portion.closed(0,n)
        elif op == ">":
            interval = portion.open(n,portion.inf)
        elif op == ">=":
            interval = portion.closedopen(n,portion.inf)
        elif op == "==":
            interval = portion.singleton(n)
        else:
            interval = portion.closedopen(1,n) | portion.open(n,portion.inf)
        return interval

    def solve(self, log=True):
        s = Solver(self)
        return s.solve(log)

    def __str__(self):
        s = ""
        for d in self.domains.values():
            s += f"{d.name}: {d}\n"
        for cf in self.choice_formulas:
            s += f"{cf}\n"
        for f in self.count_formulas:
            s += f"{f}\n"
        s += f"{self.structure}\n"
        for q in self.queries:
            s += f"query({q}).\n"
        return s
