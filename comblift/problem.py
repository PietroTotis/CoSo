
from .structure import *
from .solver import Solver
from .formulas import *
 

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
        struct_name = cof.args[0]
        cf = cof.args[1]
        op = cf.functor[1:-1]
        problog_df = cf.args[0]
        df = self.compute_dom(problog_df)
        val = cf.args[1].compute_value()
        cformula = CountingFormula(df, op, val)
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

    def add_size(self, name, s):
        self.structure.size = s

    def add_structure(self, struct):
        self.structure = struct

    def compute_dom(self, formula):
        """
        Given a formula expands and computes the corresponding domain
        Parameters
        ----------
        formula : a ProbLog predicate inter/union/not (possibly nested)
        """
        cont = self.domains[self.universe]
        if formula.functor == "inter":
            lf = self.compute_dom(formula.args[0])
            rf = self.compute_dom(formula.args[1])
            domain = lf.domain & rf.domain
        elif formula.functor == "union":
            lf = self.compute_dom(formula.args[0])
            rf = self.compute_dom(formula.args[1])
            domain = lf.domain | rf.domain
        elif formula.functor == "not":
            arg = self.compute_dom(formula.args[0])
            domain =  arg.neg().domain
        else:
            if formula.functor in self.domains:
                domain = self.domains[str(formula)]
            else:
                id = self.get_entity(formula.functor)
                if id is None:
                    raise Exception(f"Unknown constant {formula.functor}")
                domain = Domain(id, portion.singleton(id))
        df = DomainFormula(cont, formula, domain)
        return df

    def get_entity(self, e):
        if e in self.entity_map:
            return self.entity_map[e]
        else:
            return None

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
