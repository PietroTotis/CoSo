import argparse
import portion
import operator, functools
from problem import *
from formulas import *

import ply.ply.lex as lex
import ply.ply.yacc as yacc

class Lexer(object):

    def __init__(self):
        self.lexer = lex.lex(module=self)
    
    reserved = {
        'indist' : 'INDIST',
        'in' : 'IN',
        'universe' : 'UNIVERSE',
        'partitions' : 'PARTITIONS',
        'compositions' : 'COMPOSITIONS', 
    }

    # Tokens

    t_EQUALS = r'='
    t_GT = r'>'
    t_LT = r'<'
    t_COL = r':'
    t_SEMI = r';'
    # t_SLASH : r'\?'
    t_COMMA = r','
    t_LPAR  = r'\{'
    t_RPAR  = r'\}'
    t_LSPAR  = r'\['
    t_RSPAR  = r'\]'
    t_LRPAR  = r'\('
    t_RRPAR  = r'\)'
    t_COUNT = r'\#'
    t_UNION = r'\+'
    t_INTER = r'\&'
    t_DIFF  = r'!'
    t_NOT = r'\¬'
    t_ignore_COMMENT = r'%.*'
    t_ignore_WHITES = r'\ +|\t|\n'

    def t_SLASH(self, t):
        r'\|'
        t.value = "|"
        return t

    def t_LABEL(self, t):
        r'[a-z][a-zA-Z\-\_0-9]*|[0-9]+\-[0-9\-]+'
        t.type = self.reserved.get(t.value,'LABEL')    # Check for reserved words
        return t
    
    def t_NUMBER(self, t):
        r'\d+'
        try:
            t.value = int(t.value)
        except ValueError:
            print("Integer value too large %d", t.value)
            t.value = 0
        return t

    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")

    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    tokens = ['COUNT', 'COMMA', 'EQUALS', 'LT', 'GT', 'COL', 'SEMI', 'LPAR', 'RPAR', 'LSPAR', 'RSPAR', 'LRPAR', 'RRPAR', 'NUMBER', 'UNION', 'INTER', 'NOT', 'DIFF', 'LABEL', 'SLASH'] + list(reserved.values())

class Parser(object):

    tokens = Lexer.tokens

    def __init__(self, file):
        self.file = file
        self.lexer = Lexer()
        self.parser = yacc.yacc(module=self)
        self.parse_domains = True 
        self.problem = Problem()

    def parse(self):
        with open(self.file, "r") as f:
            data = f.read()
            self.parser.parse(data) # first collect domains/universe
            if len(self.problem.domains) > 0:
                self.problem.compute_universe()
                self.lexer = Lexer()
                self.parser = yacc.yacc(module=self)
                self.parse_domains = False # then structure and formula
                self.parser.parse(data)
            else:
                print("No sets found")

    def p_program(self, p):
        '''program : statement
                | statement program
        '''

    def p_statement(self, p):
        '''statement : declare_set SEMI
                | replace SEMI
                | arrangement SEMI
                | pos_constraint SEMI
                | size_constraint SEMI
        '''


    def p_entity(self, p):
        '''entity : NUMBER
                  | LABEL
        '''
        p[0] = p[1]

    def p_replace(self, p):
        '''replace : SLASH
                | SLASH SLASH
        '''
        p[0] = len(p)-1

    def p_entity_list(self, p):
        '''entity_list : entity
                        | entity COMMA entity_list
        '''
        if len(p) > 2:
            p[0] = [p[1]] + p[3]
        else:
            p[0] = [p[1]]

    def p_comp(self, p):
        '''comp : EQUALS 
            | LT
            | GT
            | GT EQUALS
            | LT EQUALS
            | DIFF EQUALS
        '''
        if len(p)>2:
            p[0] = p[1]+p[2]
        else:
            p[0] = p[1]

    def p_base_set(self, p):
        '''base_set : LABEL
               | LABEL LSPAR NUMBER RSPAR
               | UNIVERSE
        '''
        if len(p) == 2: p[0] = p[1]
        elif not self.parse_domains:
            # need more elegant way to pass number of lifted set
            p[0] = LiftedSet(self.problem.universe, SizeFormula(p[3], portion.closed(1,portion.inf)))

    def p_set(self, p):
        '''set : base_set
               | LRPAR set RRPAR
               | NOT set
               | set INTER set
               | set UNION set
        '''
        if len(p) == 2: p[0] = p[1]
        elif p[1] == '(': p[0] = p[2]
        elif p[1] == '¬': 
            if isinstance(p[1], LiftedSet):
                raise Exception("Cannot negate partition")
            else:
                p[0] = Not(p[2])
        elif p[2] == '&': 
            if isinstance(p[1], LiftedSet):
                p[1].cofs = [CountingFormula(p[3], portion.open(0,portion.inf))]
                p[0] = p[1] 
            elif isinstance(p[3], LiftedSet):
                p[3].cofs = [CountingFormula(p[1], portion.open(0,portion.inf))]
                p[0] = p[1] 
            else:
                p[0] = And(p[1],p[3])
        elif p[2] == '|': 
            if isinstance(p[1], LiftedSet):
                raise Exception("Cannot or partition")
            else:
                p[0] = Or(p[1],p[3])

    def p_declare_set(self, p):
        '''declare_set : INDIST LABEL EQUALS LPAR entity_list RPAR 
                | LABEL EQUALS LPAR entity_list RPAR 
                | INDIST LABEL EQUALS LSPAR NUMBER COL NUMBER RSPAR 
                | LABEL EQUALS LSPAR NUMBER COL NUMBER RSPAR 
        '''
        indist = p[1]=="indist"
        label = p[2] if indist else p[1]
        if ':' in p:
            lb = p[5] if indist else p[4]
            ub = p[7] if indist else p[6]
            elems = [n for n in range(lb, ub+1)]
        else:
            elems = p[5] if indist else p[4]
        ivs = self.list_to_set(elems)
        dom = functools.reduce(lambda a,b: a.union(b), ivs, portion.empty())
        dist = P.IntervalDict()
        dist[dom] = not indist
        d = DomainFormula(label, dist, self.problem.universe)
        if self.parse_domains:
            self.problem.add_domain(d)
        p[0] = d

    def p_arrangement(self, p):
        '''arrangement : LABEL IN LPAR replace set RPAR 
                       | LABEL IN LSPAR replace set RSPAR 
                       | LABEL IN PARTITIONS LRPAR set RRPAR
                       | LABEL IN COMPOSITIONS LRPAR set RRPAR 
        '''
        # not checking var consistency and set existence
        name = p[1]
        set = p[5]
        if p[3] == "partitions":
            type = "partition"   
        elif p[3] == "compositions":
            type = "composition"
        elif p[3] == '{':
            type = "subset"
        else:
            type = "sequence"
        if type == "sequence" and p[4] == 1:
            type = "permutation"
        if type == "subset" and p[4] == 2:
            type = "multisubset"
        dom = self.problem.compute_dom(set)
        s = Structure(name, type, dom)
        if self.parse_domains:
            self.problem.structure = s
        p[0] = s

    def p_sc_list(self, p):
        '''sc_list : size_constraint
                        | size_constraint COMMA entity_list
        '''
        if len(p) > 2:
            p[0] = [p[1]] + p[3]
        else:
            p[0] = [p[1]]
        print(p[0])

    def p_pos_constraint(self, p):
        '''pos_constraint : LABEL LSPAR NUMBER RSPAR EQUALS entity 
                        | LABEL LSPAR NUMBER RSPAR EQUALS set
                        | LABEL LSPAR NUMBER RSPAR EQUALS LPAR entity_list RPAR
                        | LABEL LSPAR NUMBER RSPAR EQUALS LPAR sc_list SLASH LABEL IN LABEL LSPAR NUMBER RSPAR RPAR
        '''
        arrangement = p[1]
        pos = p[3]
        if not self.parse_domains:
            if p[6] == '{':
                if p[8] == '}':
                    set = self.list_to_set(p[7])
                    df = DomainFormula("anon", set, self.problem.universe)
                else:
                    size = SizeFormula("", portion.closed(1, self.problem.universe.size()-1 ))
                    df = LiftedSet(self.problem.universe, size)
                    for sc in p[7]:
                        if isinstance(sc.formula.name, str):
                            df.size = SizeFormula("", sc.values)
                        else:
                            df.cofs.append(sc)
            else:
                df = self.problem.compute_dom(p[6])
            pf = PosFormula(arrangement, pos, df)
            self.problem.add_pos_formula(pf)
            p[0] = pf

    def p_size_constraint(self, p):
        """size_constraint : COUNT set comp NUMBER
                            | COUNT LPAR size_constraint SLASH LABEL IN LABEL RPAR comp NUMBER
        """
        p[0] = None
        if isinstance(p[2], str) and p[2] == '{':
            if not self.parse_domains:
                cf_par = p[3]
                if cf_par.formula.name == p[5]:
                    cf_par = SizeFormula(p[5], cf_par.values)
                comp = p[9]
                n = p[10]
                interval = self.problem.get_interval(comp, n)
                cf = CountingFormula(cf_par, interval)
                self.problem.add_counting_formula(cf)
                p[0] = cf
        else:
            set = p[2]
            inter = self.problem.get_interval(p[3], p[4])
            size = SizeFormula(set, inter)
            if not self.parse_domains and isinstance(set, LiftedSet):
                pos = set.size.name # yep this is ugly
                if len(set.cofs) == 0:
                    set.size = size
                else:
                    set.cofs[0].values = size
                pf = PosFormula(self.problem.structure, pos, set)
                self.problem.add_pos_formula(pf)
            else:
                if set == self.problem.structure.name:
                    if size.values.lower == 0:
                        size.values = size.values.replace(lower=1)
                    self.problem.structure.size = size
                else:
                    if not self.parse_domains:
                        df = self.problem.compute_dom(set)
                        cf = CountingFormula(df, inter)
                        self.problem.add_counting_formula(cf)
                        p[0] = cf
            


    # def p_count_constraint(self, p):
    #     '''count_constraint : COUNT set IN LABEL comp NUMBER
    #                         | COUNT LPAR size_constraint SLASH LABEL IN LABEL RPAR comp NUMBER
    #                         | COUNT LPAR count_constraint SLASH LABEL IN LABEL RPAR comp NUMBER
    #     '''
    #     if p[2] == '{':
    #         cf_par = p[3]
    #         comp = p[9]
    #         n = p[10]
    #         interval = self.problem.get_interval(comp, n)
    #         cf = CountingFormula(cf_par, interval)
    #     else:
    #         set = p[2]
    #         comp = p[5]
    #         n = p[6]
    #         interval = self.problem.get_interval(comp, n)
    #         df = self.problem.compute_dom(set)
    #         cf = CountingFormula(df, interval)
    #     self.problem.add_counting_formula(cf)
    #     p[0] = cf

    def p_error(self, p):
        if p == None:
            token = "end of file"
        else:
            token = f"{p.type}({p.value}) on line {p.lineno}"
            #ignore propery as label
            if p.type != 'PROPERTY':
                print(f"Syntax error: Unexpected {token}")

    def list_to_set(self, elems):
        entities = [self.problem.add_entity(e) for e in elems]
        entities.sort()
        ivs = []
        i = 0
        while i< len(entities):
            low = entities[i]
            while i <len(entities)-1 and entities[i]+1 == entities[i+1]: i +=1
            hi = entities[i]
            if hi - low >=1:
                ivs.append(portion.closed(low, hi))
            else:
                ivs.append(portion.singleton(low))
            i += 1
        return ivs