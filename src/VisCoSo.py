import os
import html

from yattag import Doc, indent
from util import *

HEADER = os.path.join(ROOT_DIR, "src", "VisCoSo", "header.html")

# Icons CSS

ICON_QUESTION = "fa-solid fa-question"
ICON_TIMES = "fa-solid fa-xmark"
ICON_PLUS = "fa-solid fa-plus"
ICON_MINUS = "fa-solid fa-minus"
ICON_EQUALS = "fa-solid fa-equals fa-sm"
ICON_PROBLEM = "icon fa-solid fa-pen-to-square fa-lg"
ICON_SETS = "fa-solid fa-shapes"
ICON_CONF = "fa-solid fa-cubes-stacked"
ICON_CONSTR = "fa-solid fa-ban"
ICON_INFO = "fa-solid fa-circle-info"
ICON_CALC = "fa-solid fa-calculator"

# Class names
ACCORDION = "accordion"
ACCORDIONS = "accordions"
CAPTION = "caption"
CAPTION_TEXT = "caption-text"
COLLAPSE = "collapsable"
CONTAINER = "container"
CONTENT_LINE = "content-line"
COUNT_DESC = "count-description"
COL = "column"
CONF = "configuration"
HST = "histogram"
INLINE_NUM = "inline-num"
INLINE_SOL = "inline-solution"
SUBS = "subproblems"
SOL = "solution"


class VisCoSo(object):
    def __init__(self):
        self.doc, self.tag, self.text, self.line = Doc().ttl()
        self.sol_id = 0
        self.handle = 0
        self.shatter = 0
        self.prop_colours = {}

    def add_configuration(self, type, size):
        indist = type in ["Multiset", "Set", "Partition"]
        with self.tag("table", klass=HST):
            for s in size:
                with self.tag("tr"):
                    for i in range(1, s + 1):
                        with self.tag("td", klass="slot"):
                            if not indist:
                                self.text(str(i))

    def add_content(self, problem):
        self.add_problem_repr(problem)
        self.add_solution(problem)
        has_subproblems = len(problem.subproblems) > 0
        is_shatter = len(problem.shatter_subproblems) > 0
        if has_subproblems:
            self.add_subproblems(problem.subproblems)
        if is_shatter:
            self.add_shatter(problem.shatter_subproblems)
        if not (has_subproblems or is_shatter):
            self.add_count_description(problem)

    def add_count_description(self, problem):
        with self.tag("div", klass=COLLAPSE):
            with self.tag("div", klass=CONTENT_LINE):
                pass
            with self.tag("div", klass=COUNT_DESC):
                self.icon(ICON_CALC)
                self.icon(ICON_EQUALS)
                with self.tag("p", klass=INLINE_NUM):
                    self.text(f"$${problem.count.latex()}$$")

    def add_problem(self, problem):
        with self.tag("div", klass=ACCORDION):
            with self.tag("div", klass=INLINE_SOL):
                self.icon(ICON_PROBLEM)
                with self.tag("h2", klass="problem-lab"):
                    self.text("Problem")
            with self.tag("div", klass="main-content", id="main-content"):
                self.add_content(problem)

    def add_problem_primitive(self, problem, name):
        if name == "Universe":
            cola = self.universe2cola(problem)
        elif name == "Configuration":
            cola = self.configuration2cola(problem)
        else:
            cola = self.constraints2cola(problem)

        handle = self.get_handle()
        with self.tag("div", klass=ACCORDION):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"collapse{self.handle}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass=INLINE_SOL):
                    self.icon(self.get_icon(name))
                    with self.tag("p", klass=INLINE_NUM):
                        self.text(name)
                self.label(handle)
            with self.tag("div", klass=f"{COLLAPSE} {CONF}"):
                self.visualize_primitive(name, problem)
                self.cola_collapsable(cola)

    def add_problem_repr(self, problem):
        with self.tag("div", klass=CAPTION):
            self.icon(ICON_INFO)
            with self.tag("p", klass=CAPTION_TEXT):
                self.text(problem.caption)
        with self.tag("div", klass=CONF):
            self.add_problem_primitive(problem, "Universe")
            self.add_problem_primitive(problem, "Configuration")
            self.add_problem_primitive(problem, "Constraints")

    def add_shatter(self, shatter):
        with self.tag("div", klass=COLLAPSE + " " + SUBS):
            with self.tag("div", klass=CONTENT_LINE):
                pass
            with self.tag("div", klass=ACCORDIONS):
                # with self.tag("div", klass=ACCORDION):
                with self.tag("div", klass="shatter"):
                    for pair in shatter:
                        left, right = shatter[pair]
                        op = "eq" if self.shatter == 0 else "add"
                        for r in right:
                            with self.tag("div", klass="row"):
                                self.add_subproblem(op, left)
                                self.add_subproblem("mul", r)
                                self.shatter += 1
        self.shatter = 0

    def add_subproblems(self, subproblems):
        with self.tag("div", klass=COLLAPSE + " Zero" + SUBS):
            with self.tag("div", klass=CONTENT_LINE):
                pass
            with self.tag("div", klass=ACCORDIONS):
                for i, sub in enumerate(subproblems):
                    op, subproblem = sub
                    op = "eq" if i == 0 else op
                    self.add_subproblem(op, subproblem)

    def add_subproblem(self, op, subproblem):
        handle = self.get_handle()
        with self.tag("div", klass=ACCORDION):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"collapse{self.handle}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass=INLINE_SOL):
                    self.icon(self.get_icon(op))
                    with self.tag("p", klass=INLINE_NUM):
                        self.text(str(subproblem.count.val))
                self.label(handle)
            with self.tag("div", klass=COLLAPSE):
                with self.tag("div", klass=CONTENT_LINE):
                    pass
                with self.tag("div", klass=CONTAINER):
                    self.add_content(subproblem)

    def add_solution(self, problem):
        self.sol_id += 1
        id = f"sol_{self.sol_id}"
        self.doc.input(klass="control", type="checkbox", id=id, name=id)
        with self.tag("label", ("for", id), klass=SOL):
            self.icon(ICON_EQUALS)
            with self.tag("p", klass=INLINE_NUM):
                # print(problem.caption, problem.count)
                self.text(str(problem.count.val))

    # def conf2cola(self, problem):
    #     cola_sets = ""
    #     cola_constraints = ""
    #     cola_conf = ""
    #     if problem.problem is None:
    #         cola_sets = f"universe {problem.universe.name} = {problem.universe.enumerate_elements()}; \n"
    #         for r in problem.relevant_sets:
    #             if r != problem.universe:
    #                 cola_sets += f"property {r.name} = {r.enumerate_elements()}; \n"
    #         for i, v in enumerate(problem.vars):
    #             cola_conf += f"Obj {i+1}:  {v}; \n"
    #         for c in problem.pos_constraints:
    #             cola_constraints += f"{c}; \n"
    #         for c in problem.constraints:
    #             cola_constraints += f"{c}; \n"
    #     else:
    #         p = problem.problem
    #         for d in p.domains:
    #             label = "universe" if p.domains[d] == p.universe else "property"
    #             cola_sets += f"{label} {d} = {p.domains[d].enumerate_elements()}; \n"
    #         cola_conf = str(p.configuration) + "; \n"
    #         for c in p.pos_constraints:
    #             cola_constraints += f"{c}; \n"
    #         for c in p.constraints:
    #             cola_constraints += f"{c}; \n"

    #     return cola_sets, cola_conf, cola_constraints

    def cola_collapsable(self, text):
        handle = self.get_handle()
        with self.tag("div", klass=ACCORDION):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"collapse{self.handle}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass=INLINE_SOL):
                    with self.tag("p", klass=INLINE_NUM):
                        self.text("CoLa")
                self.label(handle)
            with self.tag("div", klass=COLLAPSE):
                with self.tag("pre"):
                    self.text(text)

    def configuration2cola(self, problem):
        cola_conf = ""
        if problem.configuration is None:
            for i, v in enumerate(problem.vars):
                cola_conf += f"Obj {i+1}:  {v}; \n"
        else:
            cola_conf = str(problem.configuration)
        return cola_conf

    def constraints2cola(self, problem):
        cola_constraints = ""
        for c in problem.pos_constraints:
            cola_constraints += f"{c}; \n"
        for c in problem.constraints:
            cola_constraints += f"{c}; \n"
        return cola_constraints

    def get_handle(self):
        self.handle += 1
        return f"handle{self.handle}"

    def get_icon(self, op):
        if op == "main":
            return ICON_PROBLEM
        elif op == "mul":
            return ICON_TIMES
        elif op == "add":
            return ICON_PLUS
        elif op == "sub":
            return ICON_MINUS
        elif op == "eq":
            return ICON_EQUALS
        elif op == "Universe":
            return ICON_SETS
        elif op == "Configuration":
            return ICON_CONF
        elif op == "Constraints":
            return ICON_CONSTR
        else:
            return ICON_QUESTION

    def histogram(self, sets):
        #  relevant sets
        with self.tag("table", klass=HST):
            max = 0
            for property in sets:
                s = property.size()
                max = s if s > max else max
                if property in self.prop_colours:
                    colour_id = self.prop_colours[property]
                else:
                    colour_id = len(self.prop_colours) + 1
                    self.prop_colours[property] = colour_id
                colour = f"hc{colour_id+1}"
                with self.tag("tr"):
                    with self.tag("td", klass="ylabel"):
                        self.text(property.name)
                    for elem in property:
                        classes = f"hcell {colour}"
                        with self.tag("td", klass=classes):
                            self.text(property.universe.get_label(elem))
            if max > 4:
                with self.tag("tr"):
                    with self.tag("td", klass="ylab"):
                        self.text("")
                    for i in range(1, max + 1):
                        with self.tag("td", klass="xlab"):
                            if i % 5 == 0:
                                self.text(i)
                            else:
                                self.text("")

    def icon(self, icon):
        with self.tag("i", klass=icon):
            pass

    def label(self, handle):
        with self.tag("label", ("for", handle)):
            pass

    def visualize_primitive(self, name, problem):
        if name == "Universe":
            return self.histogram(problem.relevant_sets)
        if name == "Configuration":
            if problem.configuration is not None:
                sizes = problem.configuration.size.to_list()
            else:
                sizes = [len(problem.vars)]
            return self.add_configuration(problem.type, sizes)

    def universe2cola(self, problem):
        cola_sets = ""
        cola_sets = f"universe {problem.universe.name} = {problem.universe.enumerate_elements()}; \n"
        for r in problem.relevant_sets:
            if r != problem.universe:
                cola_sets += f"property {r.name} = {r.enumerate_elements()}; \n"
        return cola_sets

    def generate(self, log):
        h = open(HEADER, "r")
        h_content = h.read()
        h_content = html.unescape(h_content)

        self.doc.asis("<!DOCTYPE html>")
        with self.tag("html"):
            self.doc.asis(h_content)
            with self.tag("body"):
                self.add_problem(log)

        return indent(self.doc.getvalue())
