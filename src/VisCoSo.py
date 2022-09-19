import os
import html

from yattag import Doc, indent
from util import *

HEADER = os.path.join(ROOT_DIR, "src", "VisCoSo", "header.html")

ICON_QUESTION = "fa-solid fa-question"
ICON_TIMES = "fa-solid fa-xmark"
ICON_PLUS = "fa-solid fa-plus"
ICON_MINUS = "fa-solid fa-minus"
ICON_EQUALS = "fa-solid fa-equals fa-sm"
ICON_PEN = "icon fa-solid fa-pen-to-square fa-lg"
ICON_SETS = "fa-solid fa-shapes"
ICON_CONF = "fa-solid fa-cubes-stacked"
ICON_CONSTR = "fa-solid fa-ban"


class VisCoSo(object):
    def __init__(self):
        self.doc, self.tag, self.text, self.line = Doc().ttl()
        self.sol_id = 0
        self.handle = 0
        self.shatter = 0

    def add_configuration(self, problem):
        with self.tag("div", klass="caption"):
            self.icon("fa-solid fa-circle-info")
            with self.tag("p", klass="caption-text"):
                self.text(problem.caption)
        self.conf2cola(problem)

    def add_content(self, problem):
        self.add_configuration(problem)
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
        with self.tag("div", klass="collapsable"):
            with self.tag("div", klass="content-line"):
                pass
            with self.tag("div", klass="count-description"):
                self.icon("fa-solid fa-calculator")
                self.icon(ICON_EQUALS)
                with self.tag("p", klass="inline-num"):
                    self.text(f"$${problem.solution.symbolic}$$")

    def add_problem(self, problem):
        with self.tag("div", klass="accordion"):
            with self.tag("div", klass="inline-solution"):
                self.icon("icon fa-solid fa-pen-to-square fa-lg")
                with self.tag("h2", klass="problem-lab"):
                    self.text("Problem")
            with self.tag("div", klass="main-content", id="main-content"):
                self.add_content(problem)

    def add_shatter(self, shatter):
        with self.tag("div", klass="collapsable subproblems"):
            with self.tag("div", klass="content-line"):
                pass
            with self.tag("div", klass="accordions"):
                # with self.tag("div", klass="accordion"):
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
        with self.tag("div", klass="collapsable subproblems"):
            with self.tag("div", klass="content-line"):
                pass
            with self.tag("div", klass="accordions"):
                for i, sub in enumerate(subproblems):
                    op, subproblem = sub
                    op = "eq" if i == 0 else op
                    self.add_subproblem(op, subproblem)

    def add_subproblem(self, op, subproblem):
        handle = self.get_handle()
        with self.tag("div", klass="accordion"):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"collapse{self.handle}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass="inline-solution"):
                    self.icon(self.get_icon(op))
                    with self.tag("p", klass="inline-num"):
                        self.text(str(subproblem.solution.count))
                self.label(handle)
            with self.tag("div", klass="collapsable"):
                with self.tag("div", klass="content-line"):
                    pass
                with self.tag("div", klass="container"):
                    self.add_content(subproblem)

    def add_solution(self, problem):
        self.sol_id += 1
        id = f"sol_{self.sol_id}"
        self.doc.input(klass="control", type="checkbox", id=id, name=id)
        with self.tag("label", ("for", id), klass="solution"):
            self.icon(ICON_EQUALS)
            with self.tag("p", klass="inline-num"):
                self.text(str(problem.solution.count))

    def conf2cola(self, problem):
        cola_sets = ""
        cola_constraints = ""
        cola_conf = ""
        if problem.problem is None:
            cola_sets = f"universe {problem.universe.name} = {problem.universe.enumerate_elements()}; \n"
            for r in problem.relevant_sets:
                if r != problem.universe:
                    cola_sets += f"property {r.name} = {r.enumerate_elements()}; \n"
            for i, v in enumerate(problem.vars):
                cola_conf += f"Obj {i+1}:  {v}; \n"
            for c in problem.pos_constraints:
                cola_constraints += f"{c}; \n"
            for c in problem.constraints:
                cola_constraints += f"{c}; \n"
        else:
            p = problem.problem
            for d in p.domains:
                label = "universe" if p.domains[d] == p.universe else "property"
                cola_sets += f"{label} {d} = {p.domains[d].enumerate_elements()}; \n"
            cola_conf = str(p.configuration) + "; \n"
            for c in p.pos_constraints:
                cola_constraints += f"{c}; \n"
            for c in p.constraints:
                cola_constraints += f"{c}; \n"

        with self.tag("div", klass="configuration"):
            self.conf_collapsable("Sets", cola_sets)
            self.conf_collapsable("Configuration", cola_conf)
            if len(cola_constraints) > 0:
                self.conf_collapsable("Constraints", cola_constraints)

    def conf_collapsable(self, name, text):
        handle = self.get_handle()
        with self.tag("div", klass="accordion"):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"collapse{self.handle}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass="inline-solution"):
                    self.icon(self.get_icon(name))
                    with self.tag("p", klass="inline-num"):
                        self.text(name)
                self.label(handle)
            with self.tag("div", klass="collapsable"):
                with self.tag("pre"):
                    self.text(text)

    def get_handle(self):
        self.handle += 1
        return f"handle{self.handle}"

    def get_icon(self, op):
        if op == "main":
            return ICON_PEN
        elif op == "mul":
            return ICON_TIMES
        elif op == "add":
            return ICON_PLUS
        elif op == "sub":
            return ICON_MINUS
        elif op == "eq":
            return ICON_EQUALS
        elif op == "Sets":
            return ICON_SETS
        elif op == "Configuration":
            return ICON_CONF
        elif op == "Constraints":
            return ICON_CONSTR
        else:
            return ICON_QUESTION

    def icon(self, icon):
        with self.tag("i", klass=icon):
            pass

    def label(self, handle):
        with self.tag("label", ("for", handle)):
            pass

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
