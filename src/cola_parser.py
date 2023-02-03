import portion as Int
import ply.lex as lex
import ply.yacc as yacc

from .problem import Problem
from .venn import Venn
from .configuration import *
from .level_1 import *
from .level_2 import *


class Lexer(object):
    def __init__(self):
        self.lexer = lex.lex(module=self)

    reserved = {
        "labelled": "LAB",
        "property": "PROP",
        "in": "IN",
        "universe": "UNIVERSE",
        "part": "PART",
        "repeated": "REPEAT",
    }

    # Tokens

    t_EQUALS = r"="
    t_GT = r">"
    t_LT = r"<"
    # t_COL = r":"
    t_SEMI = r";"
    # t_SLASH : r'\?'
    t_COMMA = r","
    t_LPAR = r"\{"
    t_RPAR = r"\}"
    t_LSPAR = r"\["
    t_RSPAR = r"\]"
    t_LRPAR = r"\("
    t_RRPAR = r"\)"
    t_COUNT = r"\#"
    t_UNION = r"\+"
    t_INTER = r"\&"
    t_DIFF = r"!"
    t_NOT = r"\¬"
    t_ignore_COMMENT = r"%.*"
    t_ignore_WHITES = r"\ +|\t|\n"

    def t_LABEL(self, t):
        r"[a-z][a-zA-Z\-\_0-9]*|[0-9]+\-[0-9\-]+"
        t.type = self.reserved.get(t.value, "LABEL")  # Check for reserved words
        return t

    def t_NUMBER(self, t):
        r"\d+"
        try:
            t.value = int(t.value)
        except ValueError:
            print("Integer value too large %d", t.value)
            t.value = 0
        return t

    def t_newline(self, t):
        r"\n+"
        t.lexer.lineno += t.value.count("\n")

    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    tokens = [
        "COUNT",
        "COMMA",
        "EQUALS",
        "LT",
        "GT",
        # "COL",
        "SEMI",
        "LPAR",
        "RPAR",
        "LSPAR",
        "RSPAR",
        "LRPAR",
        "RRPAR",
        "NUMBER",
        "UNION",
        "INTER",
        "NOT",
        "DIFF",
        "LABEL",
    ] + list(reserved.values())


class Parser(object):

    tokens = Lexer.tokens

    def __init__(self, file=None, cola=None):
        self.file = file
        self.cola = cola
        self.lexer = Lexer()
        self.parser = yacc.yacc(module=self)
        self.sizes = {}
        self.problem = Problem()
        self.venn = Venn()

    def parse(self):
        if self.file is not None:
            with open(self.file, "r") as f:
                cola = f.read()
        else:
            cola = self.cola
        self.parser.parse(cola)
        if self.problem.configuration is None:
            raise EmptyException("No configuration specified")

    def p_program(self, p):
        """program : statement
        | statement program
        """

    def p_statement(self, p):
        """statement : declare_set SEMI
        | arrangement SEMI
        | pos_constraint SEMI
        | count_constraint SEMI
        """

    def p_entity_list(self, p):
        """entity_list : LABEL
        | LABEL COMMA entity_list
        """
        if len(p) > 2:
            p[0] = [p[1]] + p[3]
        else:
            p[0] = [p[1]]

    def p_comp(self, p):
        """comp : EQUALS
        | LT
        | GT
        | GT EQUALS
        | LT EQUALS
        | DIFF EQUALS
        """
        if len(p) > 2:
            p[0] = p[1] + p[2]
        else:
            p[0] = p[1]

    def p_base_set(self, p):
        """base_set : LABEL
        | LABEL LSPAR NUMBER RSPAR
        | UNIVERSE
        """
        if len(p) == 2:
            p[0] = p[1]
        else:
            # need more elegant way to pass number of lifted set
            p[0] = LiftedSet(self.problem.universe, CSize(p[3], Int.closed(1, Int.inf)))

    def p_set_atom(self, p):
        """set_atom : base_set
        | LRPAR set RRPAR
        """
        if p[1] == "(":
            p[0] = p[2]
        else:
            p[0] = p[1]

    def p_set_not(self, p):
        """set_not : NOT set_atom
        | set_atom
        """
        p[0] = p[1]

    def p_set_inter(self, p):
        """set_inter : set_not
        | set_not INTER set_atom
        """
        if len(p) > 2:
            if isinstance(p[1], LiftedSet):
                d = self.problem.compute_dom(p[3])
                p[1].ccs = [CCounting(d, Int.open(0, Int.inf))]
                p[0] = p[1]
            elif isinstance(p[3], LiftedSet):
                d = self.problem.compute_dom(p[1])
                p[3].ccs = [CCounting(d, Int.open(0, Int.inf))]
                p[0] = p[1]
            else:
                p[0] = And(p[1], p[3])
        else:
            p[0] = p[1]

    def p_set_union(self, p):
        """set_union : set_inter
        | set UNION set_inter
        """
        if len(p) > 2:
            if isinstance(p[1], LiftedSet):
                raise Exception("Cannot or partition")
            else:
                p[0] = Or(p[1], p[3])
        else:
            p[0] = p[1]

    def p_set(self, p):
        """set : set_union
        | PART
        """

        if len(p) == 2:
            p[0] = p[1]
        else:
            if isinstance(p[2], LiftedSet):
                raise Exception("Cannot negate partition")
            else:
                p[0] = Not(p[2])

    def p_declare_set(self, p):
        """declare_set : PROP LABEL EQUALS LPAR entity_list RPAR
        | UNIVERSE LABEL EQUALS LPAR entity_list RPAR
        | PROP UNIVERSE EQUALS LPAR entity_list RPAR
        | LAB PROP LABEL EQUALS LPAR entity_list RPAR
        | PROP LABEL
        | LAB PROP LABEL
        """
        enumerated = "=" in p
        implicit_labelled = p[1] == "labelled"
        explicit_univ = p[1] == "universe"
        distinguishable = enumerated or implicit_labelled
        label = p[3] if implicit_labelled else p[2]
        if "=" not in p:
            self.venn.add_base_set(label, distinguishable)
        else:
            k = 1 if implicit_labelled else 0
            elems = p[5 + k]
            implicit_univ = enumerated and self.problem.universe.empty()
            univ = implicit_univ or explicit_univ
            elems = list2interval(self.problem, elems, univ)
            if explicit_univ:
                d = Universe(label, elems, name=label)
                self.problem.universe = d
            else:
                d = SetFormula(label, elems, self.problem.universe, name=label)
            self.problem.add_domain(d)
            p[0] = d

    def p_arrangement(self, p):
        """arrangement : LABEL IN LPAR REPEAT set RPAR
        | LABEL IN LPAR set RPAR
        | LABEL IN LSPAR REPEAT set RSPAR
        | LABEL IN LSPAR set RSPAR
        | LABEL IN LPAR LPAR set RPAR RPAR
        | LABEL IN LSPAR LPAR set RPAR RSPAR
        """
        # not checking var consistency and set existence
        name = p[1]
        if p[3] == "[":
            if p[4] == "repeated":
                type = "sequence"
            elif p[4] == "{":
                type = "composition"
            else:
                type = "permutation"
        else:  # p[3] == "{"
            if p[4] == "repeated":
                type = "multisubset"
            elif p[4] == "{":
                type = "partition"
            else:
                type = "subset"
        set = p[4] if type in ["subset", "permutation"] else p[5]
        if len(self.venn.base_sets) > 0:
            for set in self.sizes:
                self.venn.add_set(set, self.sizes[set])
            self.problem = self.venn.update_domains(self.problem)
        dom = self.problem.compute_dom(set)
        s = Configuration(name, type, dom)
        self.problem.configuration = s
        if self.problem.universe.empty():
            self.problem.compute_universe()
        if len(self.problem.domains) > 0:
            names = list(self.problem.domains)
            for i, d1 in enumerate(names):
                for d2 in names[i:]:
                    dom1 = self.problem.domains[d1]
                    dom2 = self.problem.domains[d2]
                    if dom1.elements.domain() in dom2.elements.domain():
                        dom1.elements.combine(
                            dom2.elements, how=is_distinguishable
                        )  # update indistinguishable
                    if dom2.elements.domain() in dom1.elements.domain():
                        self.problem.domains[d1] = dom1 | dom2
        else:
            raise EmptyException("No sets found")
        p[0] = s

    def p_pos_constraint(self, p):
        """pos_constraint : LABEL LSPAR NUMBER RSPAR EQUALS set
        | LABEL LSPAR NUMBER RSPAR IN set
        | LABEL LSPAR NUMBER RSPAR EQUALS LPAR entity_list RPAR
        """
        arrangement = p[1]
        pos = p[3]
        if p[6] == "{":
            if p[8] == "}":
                set = list2interval(self.problem, p[7], False)
                df = SetFormula("anon", set, self.problem.universe)
            else:
                size = CSize("", Int.closed(1, self.problem.universe.size() - 1))
                df = LiftedSet(self.problem.universe, size)
                for sc in p[7]:
                    if isinstance(sc.formula.name, str):
                        df.size = CSize("", sc.values)
                    else:
                        df.ccs.append(sc)
        else:
            df = self.problem.compute_dom(p[6])
        pf = CPosition(arrangement, pos, df)
        self.problem.add_pos_formula(pf)
        p[0] = pf

    def p_count_constraint(self, p):
        """count_constraint : COUNT set comp NUMBER
        | COUNT LRPAR count_constraint RRPAR comp NUMBER
        """
        if isinstance(p[2], str) and p[2] == "(":
            # level 2 counting constraint
            cf_par = p[3]
            if cf_par.formula.name == self.problem.universe.name:
                cf_par = CSize(p[5], cf_par.values)
            comp = p[5]
            n = p[6]
            interval = self.problem.get_interval(comp, n)
            cf = CCounting(cf_par, interval)
            self.problem.add_counting_formula(cf)
            p[0] = cf
        else:  # level 1 counting constraint
            set = p[2]
            inter = self.problem.get_interval(p[3], p[4])
            size = CSize(set, inter)
            if isinstance(set, LiftedSet):
                # create a lifted set for each position and later conjoin with
                # variable
                pos = set.size.name  # yep this is ugly
                if len(set.ccs) == 0:
                    set.size = size
                else:
                    set.ccs[0].values = inter
                pf = CPosition(self.problem.configuration, pos, set)
                self.problem.add_pos_formula(pf)
            else:
                conf = self.problem.configuration
                if conf is None:
                    # parsing sizes of sets before config declaration
                    self.venn.add_set(set, p[4])
                elif conf is not None and set == conf.name:
                    if size.values.lower == 0:
                        size.values = size.values.replace(lower=1)
                    self.problem.configuration.size = size
                else:
                    df = self.problem.compute_dom(set)
                    cf = CCounting(df, inter)
                    self.problem.add_counting_formula(cf)
                    p[0] = cf

    def p_error(self, p):
        if p == None:
            token = "end of file"
        else:
            token = f"{p.type}({p.value}) on line {p.lineno}"
            # ignore propery as label
            if p.type != "PROPERTY":
                print(f"Syntax error: Unexpected {token}")
