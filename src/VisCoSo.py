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
        self.line("pre", self.conf2cola(problem))

    def add_content(self, problem):
        self.add_configuration(problem)
        self.add_solution(problem)
        has_subproblems = len(problem.subproblems) > 0
        is_shatter = len(problem.shatter_subproblems) > 0
        if has_subproblems:
            self.add_subproblems(problem.subproblems)
        if is_shatter:
            self.add_shatter(problem.shatter_subproblems, problem.solution.count)
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

    def add_shatter(self, shatter, solution):
        handle = self.get_handle()
        with self.tag("div", klass="accordion"):
            self.doc.input(
                klass="control",
                type="checkbox",
                id=handle,
                name=f"sol_{self.sol_id}",
            )
            with self.tag("div", klass="handle"):
                with self.tag("div", klass="inline-solution"):
                    self.icon(self.get_icon("eq"))
                    with self.tag("p", klass="inline-num"):
                        self.text(str(solution))
                self.label(handle)
            with self.tag("div", klass="collapsable"):
                with self.tag("div", klass="content-line"):
                    pass
                for pair in shatter:
                    left, right = shatter[pair]
                    op = "eq" if self.shatter == 0 else "add"
                    for r in right:
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
        cola = f"universe u = {problem.universe};\n"
        for i, v in enumerate(problem.vars):
            cola += f"Obj {i+1}:  {v};\n"
        for c in problem.pos_constraints:
            cola += f"{c};\n"
        for c in problem.constraints:
            cola += f"{c};\n"
        return cola

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
