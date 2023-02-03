import itertools
from .VisCoSo import VisCoSo
from .level_1 import Multiset
from .level_2 import LiftedSet
from .configuration import CCounting
from .util import *


class ActionLog(object):
    def __init__(self, description):
        self.description = description
        self.details = []

    def detail(self, msg):
        self.details.append(msg)


class ProblemLog(object):
    def __init__(
        self,
        vars=[],
        config=None,
        pos_constraints=[],
        constraints=[],
        universe=None,
        level=0,
        id="1",
        debug=True,
        caption="Problem",
        configuration=None,
    ):
        self.actions = []
        self.caption = caption
        self.configuration = configuration
        self.constraints = constraints
        self.debug = debug
        self.id = id
        self.level = level
        self.pos_constraints = pos_constraints
        self.relevant_sets = [universe]
        self.shatter_cases = {}
        self.shatter_subproblems = {}
        self.count = None
        self.space = "  "  # "\t"
        self.subproblems = []
        self.config = config
        self.vars = [v.copy() for v in vars]
        self.universe = universe
        self.vis = VisCoSo()
        self.indent = f"{self.space}" * self.level
        for c in self.pos_constraints:
            self.add_relevant_set(c.formula)
        for c in self.constraints:
            self.add_relevant_set(c.formula)
        if debug:
            print(self.configuration2text())

    def __str__(self):
        msg = self.configuration2text() + "\n"
        msg += self.actions2text()
        msg += self.subproblems2text()
        msg += self.shatter2text()
        msg += self.indent + "========\n"
        msg += f"{self.indent}({self.id}) Solution: {self.count}\n"
        return msg

    def copy(self):
        copy = ProblemLog(
            vars=self.vars,
            config=self.config,
            pos_constraints=self.pos_constraints,
            constraints=self.constraints,
            universe=self.universe,
            level=self.level,
            id=self.id,
            debug=self.debug,
            caption=self.caption,
            configuration=self.configuration,
        )
        return copy

    def action(self, msg):
        al = ActionLog(msg)
        self.actions.append(al)
        if self.debug:
            print(f"{self.indent}- {msg}")
        return al

    def add_relevant_set(self, domain):
        if isinstance(domain, Multiset):
            self.relevant_sets = self.relevant_cases(domain, self.relevant_sets)
        elif isinstance(domain, CCounting):  # level 2 constraint
            self.relevant_sets = self.relevant_cases(domain.formula, self.relevant_sets)
        elif isinstance(domain, LiftedSet):  # level 2 pos constraint
            for c in domain.ccs:
                self.add_relevant_set(c.formula)
        else:  # ignore CSize
            pass

    def add_subproblem(self, op, sub_log):
        """Append subproblem

        Args:
            op (str): descriptor for the operation aggregating the subproblems, i.e. add/mul/...
            sub_log (ProblemLog): _description_
        """
        self.subproblems.append((op, sub_log))
        if self.debug:
            print(f"{sub_log.indent}========")
            print(f"{sub_log.indent}({sub_log.id}) Solution: {sub_log.count}")

    def add_split_left(self, subproblem, shatter_id):
        self.shatter_subproblems[shatter_id] = (subproblem, [])
        if self.debug and subproblem.solution is not None:
            print(f"{subproblem.indent}========")
            print(f"{subproblem.indent}({subproblem.id}) Solution: {subproblem.count}")

    def add_split_right(self, subproblem, shatter_id):
        _, right = self.shatter_subproblems[shatter_id]
        right.append(subproblem)
        if self.debug:
            print(f"{subproblem.indent}========")
            print(f"{subproblem.indent}({subproblem.id}) Solution: {subproblem.count}")

    def detail(self, action, msg):
        if isinstance(msg, str):
            action.detail(msg)
            if self.debug:
                print(f"{self.indent}{self.space} {msg}")
        elif isinstance(msg, list):
            for elem in msg:
                action.detail(str(elem))
                if self.debug:
                    print(f"{self.indent}{self.space} {elem}")

    def description(self, msg):
        self.caption = msg

    def propagation(self, constraint):
        al = self.action("Propagate")
        al.detail(str(constraint))
        if self.debug:
            print(f"{self.indent}- Propagating {constraint}")
        return al

    def shatter_case(self, action, n, split, rest):
        self.shatter_cases[n] = (split, rest)
        if self.debug:
            print(f"{self.indent}C{n}) {split} // {rest}")

    def warning(self, msg):
        al = ActionLog("warning")
        al.detail(msg)

    #######################
    ### Text formatting ###
    #######################

    def configuration2text(self):
        text = (
            self.indent + f"### ({self.id}) {self.caption} ({self.config.type}) ###\n"
        )
        if len(self.vars) > 0:
            for i, v in enumerate(self.vars):
                text += self.indent + f"Obj {i+1}:  {v}\n"
        if len(self.pos_constraints) > 0:
            text += self.indent + "Positional constraints:\n"
            for c in self.pos_constraints:
                text += self.indent + f"{self.space}{c}\n"
        if len(self.constraints) > 0:
            text += self.indent + "Counting constraints:\n"
            for c in self.constraints:
                text += self.indent + f"{self.space}{c}\n"
        text += self.indent + "%%%%%%%"
        return text

    def action2text(self, action):
        text = self.indent + f"- {action.description}\n"
        if len(action.details) > 0:
            text += self.indent + "--------\n"
            for detail in action.details:
                text += self.indent + self.space + detail + "\n"
        return text

    def actions2text(self):
        text = ""
        for action in self.actions:
            text += self.action2text(action)
        if len(self.actions) > 0:
            text += self.indent + "~~~~~~~~\n"
        return text

    def relevant_cases_intersection(self, mset, rest_classes):
        """
        Computes recursively the intersections of relevant cases/domains
        """
        base = [mset, mset.neg()]
        if len(rest_classes) == 0:
            return base
        else:
            combinations = [base]
            first = rest_classes[0]
            if first == self.universe:
                return self.relevant_cases_intersection(mset, rest_classes[1:])
            for dom in rest_classes[1:]:
                if not dom == self.universe:
                    comb = [dom, dom.neg()]
                    combinations.append(comb)
            combinations = list(itertools.product(*combinations))
            cases = []
            for c in combinations:
                dom_base = c[0]
                if len(c) > 1:
                    for dom in c[1:]:
                        dom_base = dom_base & dom
                if dom_base.size() > 0:
                    cases.append(dom_base)
            return cases

    def relevant_cases(self, mset, rest_classes):
        if len(rest_classes) == 0:
            res = [mset, mset.neg()]
        else:
            res = self.relevant_cases_intersection(mset, rest_classes)
        # print(res)
        return res

    def subproblems2text(self):
        text = ""
        for op, subproblem in self.subproblems:
            if op == "add":
                text += f"{self.indent}Summing:\n"
            elif op == "mul":
                text += f"{self.indent}Multiplying:\n"
            elif op == "sub":
                text += f"{self.indent}Subtracting:\n"
            else:
                text += f"{self.indent}Subproblem:\n"
            text += str(subproblem)
        return text

    def shatter2text(self):
        text = ""
        if len(self.shatter_subproblems) > 0:
            text = f"{self.indent}Plus\n"
            for id in self.shatter_subproblems:
                left, right = self.shatter_subproblems[id]
                text += str(left)
                if len(right) > 0:  # right might be missing if left=0
                    text += f"{self.indent}Times\n"
                    for rest_subproblem in right:
                        text += str(rest_subproblem)
        return text

    def to_viscoso(self):
        return self.vis.generate(self)

    def to_viscoso_widget(self):
        return self.vis.generate_widget(self)
