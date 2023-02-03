import os
import html
import random
import portion as P
from yattag import Doc, indent

from .level_1 import Multiset
from .level_2 import LiftedSet
from .configuration import CSize
from .util import *

HEADER = os.path.join(ROOT_DIR, "src", "VisCoSo", "header.html")

# Icons CSS

ICON_QUESTION = "fa-solid fa-question"
ICON_TIMES = "fa-solid fa-times"
ICON_PLUS = "fa-solid fa-plus"
ICON_MINUS = "fa-solid fa-minus"
ICON_EQUALS = "fa-solid fa-equals fa-sm"
ICON_PROBLEM = "fa-solid fa-file-alt"
ICON_SETS = "fa-solid fa-shapes"
ICON_CONF = "fa-solid fa-sitemap"
ICON_CONSTR = "fa-solid fa-ban"
ICON_INFO = "fa-solid fa-info-circle"
ICON_CALC = "fa-solid fa-calculator"
ICON_ARROWL = "fa-solid fa-arrow-left"
ICON_ARROWR = "fa-solid fa-arrow-right"
ICON_NOPE = "fa-solid fa-x"
ICON_GROUP = "fa-solid fa-object-group"

# Class names
ACCORDION = "accordion"
ACCORDIONS = "accordions"
CAPTION = "caption"
CAPTION_TEXT = "caption-text"
COLLAPSE = "collapsable"
CONTAINER = "vs-container"
CONTENT_LINE = "content-line"
COUNT_DESC = "count-description"
COL = "column"
CONF = "configuration"
HST = "histogram"
INLINE_NUM = "inline-num"
INLINE_SOL = "inline-solution"
LSET = "lifted-set"
SUBS = "subproblems"
SOL = "solution"

# CSS styles
COLOURS = {
    "red": "#bf616a",
    "orange": "#d08770",
    "yellow": "#ebcb8b",
    "green": "#a3be8c",
    "violet": "#b48ead",
    "blue": "#5e81ac",
    "lightblue": "#81a1c1",
    "cyan": "#88c0d0",
    "lightgreen": "#8fbcbb",
    "none": "#4c566a",
}


def fill(colour):
    return f"background-color: {colour}"


def border_col(colour):
    return f"border-color: {colour}"


class VisCoSo(object):
    def __init__(self):
        self.doc, self.tag, self.text, self.line = Doc().ttl()
        self.sol_id = 0
        self.handle = 0
        self.shatter = 0
        self.prop_colours = {}

    def add_configuration(self, config, size, vars):
        if config.lvl1():
            self.add_lvl1_configuration(config, size, vars)
        if config.lvl2():
            self.add_lvl2_configuration(config, size, vars)

    def add_lvl1_configuration(self, config, size, vars):
        with self.tag("table", klass=HST):
            for s in size:
                with self.tag("tr"):
                    for i in range(0, s):
                        colour = "transparent"
                        if len(vars) > 0 and vars[i] in self.prop_colours:
                            colour = self.get_colour(vars[i], "transparent")
                        with self.tag("td", klass="slot", style=f"{fill(colour)}"):
                            if config.labelled():
                                self.text(str(i + 1))

    def add_lvl2_configuration(self, config, size, vars):
        with self.tag("table", klass=HST):
            for s in size:
                with self.tag("tr slot"):
                    for i in range(0, s):
                        if len(vars) > 0:
                            self.add_lifted_set(vars[i], labelled=config.labelled)
                        else:
                            ub = config.domain.size() - s + 1
                            self.add_lifted_set(
                                max=ub, labelled=config.labelled, id=i + 1
                            )
            # if config.labelled():
            #     with self.tag("tr slot"):
            #         for i in range(0, size[-1]):
            #             with self.tag("td"):
            #                 self.text(str(i + 1))

    def add_constraints(self, constraints):
        with self.tag("table", klass=HST):
            for c in constraints:
                with self.tag("tr"):
                    with self.tag("td"):
                        self.icon(ICON_CONF)
                    with self.tag("td"):
                        self.icon(ICON_ARROWR)
                    self.add_constraint(c)

    def add_constraint(self, c):
        if isinstance(c.formula, Multiset):
            self.add_lvl1_constraint(c)
        else:
            self.add_lvl2_constraint(c)

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

    def add_interval(self, inter, colour, box_type="constr-cell"):
        vals = list([i for i in P.iterate(inter, step=1)])
        # print(inter, vals, vals[-1])
        filler = fill(colour) if box_type == "constr-cell" else ""
        if vals == [0]:
            with self.tag("td", klass=f"{box_type}", style=f"{border_col(colour)}"):
                self.icon(ICON_NOPE)
        for s in range(0, vals[0]):
            with self.tag(
                "td", klass=f"{box_type}", style=f"{filler}; {border_col(colour)}"
            ):
                pass
        for s in range(vals[0], vals[-1]):
            if s in vals:
                if s != 0:
                    with self.tag(
                        "td", klass=f"{box_type}-maybe", style=f"{border_col(colour)}"
                    ):
                        pass
            else:
                with self.tag("td", klass=f"{box_type}", style=f"{border_col(colour)}"):
                    self.icon(ICON_NOPE)

    def add_lifted_set(self, lset=None, labelled=False, max=None, id=0):
        if labelled:
            name = f"Part {id}" if lset is None else lset.name
        else:
            name = f"Part" if lset is None else lset.name
        sizes = P.closed(1, max) if lset is None else lset.size.values
        ccs = [] if lset is None else lset.ccs
        histogram = {} if lset is None else lset.histogram
        max = max if lset is None else lset.size.values.upper

        with self.tag("table", klass=LSET):
            with self.tag("tr"):
                with self.tag("td", colspan=f"{max+1}", style="text-align: center"):
                    self.icon(ICON_GROUP)
                    self.text(name)
            with self.tag("tr"):
                with self.tag("td"):
                    self.text("Size:")
                self.add_interval(sizes, COLOURS["none"])
            for cc in ccs:
                if cc.formula not in histogram:
                    with self.tag("tr"):
                        with self.tag("td"):
                            self.text(str(cc.formula))
                            self.add_constraint(cc)
            for prop in histogram.keys():
                val = histogram[prop]
                if val > -1:
                    with self.tag("tr"):
                        with self.tag("td"):
                            self.text(str(prop))
                            self.add_interval(P.singleton(val), self.get_colour(prop))

    def add_lvl1_constraint(self, c):
        if isinstance(c, CSize):
            colour = COLOURS["none"]
        else:
            colour = self.get_colour(c.formula)
        self.add_interval(c.values, colour)

    def add_lvl2_constraint(self, c):
        self.add_interval(c.values, COLOURS["none"], box_type="slot")
        with self.tag("td"):
            self.icon(ICON_ARROWL)
        self.add_lvl1_constraint(c.formula)

    def add_pos_constraints(self, pcs):
        with self.tag("table", klass=HST):
            for pc in pcs:
                prop = pc.formula
                colour = self.get_colour(prop)
                with self.tag("tr"):
                    with self.tag("td", klass="slot"):
                        self.text(str(pc.pos))
                    with self.tag("td"):
                        self.icon(ICON_ARROWL)
                    with self.tag("td", klass=f"hcell", style=f"{fill(colour)}"):
                        pass
                    with self.tag("td", klass="ylabel"):
                        self.text(str(prop))

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
            if len(problem.pos_constraints) > 0 or len(problem.constraints) > 0:
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

    def get_colour(self, property, default=None):
        if property.name in COLOURS:
            return COLOURS[property.name]
        else:
            if property in self.prop_colours:
                return self.prop_colours[property]
            elif default is not None:
                return default
            else:
                available = set(COLOURS.values()) - set(self.prop_colours.values())
                if available:
                    colour = available.pop()
                else:
                    colour = hex(random.randrange(0, 2**24))
                self.prop_colours[property] = colour
                return colour

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
                colour = self.get_colour(property)
                with self.tag("tr"):
                    with self.tag("td", klass="ylabel"):
                        self.text(property.name)
                    for elem in property:
                        with self.tag("td", klass="hcell", style=f"{fill(colour)}"):
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
        with self.tag("i", klass=f"vsc-icon {icon}"):
            pass

    def label(self, handle):
        with self.tag("label", ("for", handle), klass="vsc-lab"):
            pass

    def reserve_colours(self, log):
        for rvset in log.relevant_sets:
            if rvset.name in COLOURS:
                self.prop_colours[rvset] = COLOURS[rvset.name]

    def visualize_primitive(self, name, problem):
        if name == "Universe":
            return self.histogram(problem.relevant_sets)
        elif name == "Configuration":
            if problem.configuration is not None:
                sizes = problem.configuration.size.to_list()
            else:
                sizes = [len(problem.vars)]
            return self.add_configuration(problem.config, sizes, problem.vars)
        else:  # Constraint
            self.add_pos_constraints(problem.pos_constraints)
            self.add_constraints(problem.constraints)

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

        self.reserve_colours(log)

        self.doc.asis("<!DOCTYPE html>")
        with self.tag("html"):
            self.doc.asis(h_content)
            with self.tag("body"):
                self.add_problem(log)

        return indent(self.doc.getvalue())

    def generate_widget(self, log):
        self.reserve_colours(log)
        with self.tag(
            "script",
            src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML",
            type="text/javascript",
        ):
            pass
        with self.tag("script", type="text/x-mathjax-config"):
            self.text("""MathJax.Hub.Config({TeX: {extensions: ["action.js"] }});""")
        with self.tag("div"):
            with self.tag("div"):
                self.text("$$\\require{action}$$")  # make sure texttips are active
            self.add_problem(log)
        return indent(self.doc.getvalue())
