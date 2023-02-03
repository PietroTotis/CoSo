from portion.dict import IntervalDict

from .level_1 import SetFormula
from .util import *

class Union(object):
    def __init__(self, formula, size=None, subsets=[]):
        self.formula = formula
        self.size = size
        self.subsets = subsets

    def __repr__(self):
        descr = f"{self.formula}: {self.size}\n"
        for s in self.subsets:
            descr += f"\t{s}\n"
        return descr


class Venn(object):
    """
    A class to represent a Venn diagram: the universe is divided into subsets
    associated to the number of objects belonging to the subsets.

    """

    def __init__(self):
        self.base_sets = []
        self.base_indist = {}
        self.indist = {}
        self.sets = {}
        self.parts = {}

    def __repr__(self):
        descr = ""
        for set in self.sets:
            descr += f"#({self.area2formula(set)}) = {self.sets[set]}\n"
        return descr

    def add_base_set(self, set, dist):
        """
        Add a set which is a label and not a formula

        Args:
            set (str): label
            dist (bool): the elements are labelled
        """
        self.base_sets.append(set)
        self.base_indist[set] = not dist

    def add_set(self, set, size):
        dnfied = dnfy(set)
        if isinstance(dnfied, Or):
            subsets = [s for s in flatten(Or, dnfied)]
        else:
            subsets = []
        u = Union(dnfied, size, subsets)
        self.sets[set] = u

    def subsets(self):
        """
        For all base sets collect subsets from declared intersections
        """
        sets = [self.formula2area(b) for b in self.base_sets]
        for set in sets:
            subs = self.get_subsets(set)
            f = self.sets.get(self.area2formula(set))
            if f is not None:
                s = f.size
            else:
                s = None
            name = self.area2formula(set)
            if len(subs) > 0:
                if name not in self.sets:
                    self.sets[name] = Union(name, s, subs)
                else:
                    self.sets[name].subsets += subs
            else:
                self.sets[name] = Union(name, s, [])
                # self.parts.append(set)

    def get_subsets(self, set1):
        """
        Check which of the sets with declared size are subsets of set 1
        """
        declared = [self.formula2area(k) for k in self.sets.keys()]
        subsets = set()
        for set2 in declared:
            if set1 != set2:
                if self.is_subset(set2, set1):
                    subsets2 = self.get_subsets(set2)
                    if len(subsets2) > 0:
                        subsets.union(subsets2)
                    else:
                        subsets.add(set2)
        return subsets

    def build_parts(self):
        univ_area = (0,) * len(self.base_sets)
        parts = self.get_parts(univ_area)
        for p in parts:
            f = self.area2formula(p)
            if f in self.sets:
                self.parts[p] = self.sets[f].size
        return parts

    def is_subset(self, set1, set2):
        """
        Checks if set1 is contained in set2 by comparing codes

        Args:
            set1 (tuple(int)): code of set 1
            set2 (tuple(int)): code of set 2

        Returns:
            bool: set1 in set2
        """
        subset = True
        s1 = self.pad(set1)
        s2 = self.pad(set2)
        for i in range(0, len(self.base_sets)):
            if s1[i] != 0 and s1[i] == -s2[i]:
                subset = False
            if s1[i] == 0 and s2[i] == 1:
                subset = False
        return subset

    def infer_union_size(self, u):
        """
        Inference on the sizes of the sets

        Args:
            u (Union): a set

        Raises:
            Exception: we might find inconsistent sizes
            Exception: we might miss some sizes
        """
        f = self.formula2area(u.formula)
        parts = [tuple(p) for p in self.get_parts(f)]
        known = [p for p in parts if p in self.parts]
        ksize = sum([self.parts[part] for part in known])
        unknown = [p for p in parts if p not in self.parts]
        if len(unknown) == 0:
            if u.size is None:
                u.size = ksize
            elif u.size != ksize:
                raise Exception(
                    f"Inconsistent size of {u.formula} (expected {u.size}, got {ksize})"
                )
        elif len(unknown) == 1:
            if u.size is None:
                u.size = ksize
                self.parts[unknown[0]] = 0
            else:
                self.parts[unknown[0]] = u.size - ksize
        else:
            if u.size is None:
                if len(known) == 0:
                    raise Exception(
                        f"I don't know neither the size of {u.formula} nor any of its subsets"
                    )
                else:
                    u.size = ksize
            else:
                if len(known) == 0:
                    self.parts[f] = u.size
                # missing the case where both #A and #B are unknown but #A u B is known

    def infer_indistinguishability(self):
        for p in self.parts:
            indist = self.parts[p] > 1
            sets = [self.base_sets[i] for i, s in enumerate(p) if s == 1]
            for set in sets:
                indist = indist and self.base_indist[set]
            self.indist[p] = indist

    def infer(self):
        self.subsets()
        self.build_parts()
        for set in self.sets:
            self.infer_union_size(self.sets[set])
        self.infer_indistinguishability()

    def area2formula(self, area):
        """
        Convert an area code int

        Args:
            area (tuple(-1/1/0)): whether a base set is excluded/included/neither

        Returns:
            And: intersection of base sets or their negation according to the code
        """
        bsets = []
        for i, n in enumerate(area):
            if n == 1:
                bsets.append(self.base_sets[i])
            elif n == -1:
                bsets.append(Not(self.base_sets[i]))
            else:
                pass
        return nest(And, bsets)

    def formula2area(self, formula):
        part = [0] * len(self.base_sets)
        sets = flatten(And, formula)
        for s in sets:
            negated = isinstance(s, Not)
            if negated:
                i = self.base_sets.index(s.child)
                part[i] = -1
            else:
                i = self.base_sets.index(s)
                part[i] = 1
        return tuple(part)

    def pad(self, area):
        """
        Add zeros to a set description created before adding all base sets

        Args:
            area (tuple(int)): the partial encoding
        Returns:
            tuple(int): same encoding with additional zeros
        """
        padded = []
        for id in area:
            padded.append(id)
        pad = len(self.base_sets) - len(area)
        padded += [0] * pad
        return tuple(padded)

    def get_parts(self, area):
        """
        Given a region code of the Venn diagram returns all the partition
        codes that the area covers

        Args:
            area (tuple(int)): code of the set

        Returns:
            [tuple(int)]: list of parts in the set
        """
        if len(area) == 0:
            return [()]
        else:
            try:
                z = area.index(0)
                p1 = [area[:z] + (1,) + r for r in self.get_parts(area[z + 1 :])]
                p2 = [area[:z] + (-1,) + r for r in self.get_parts(area[z + 1 :])]
                return p1 + p2
            except ValueError:
                return [area]

    def update_domains(self, problem):
        """
        Convert a Venn diagram into domains with the correct intersections

        Args:
            problem (Problem): problem corresponding to the Venn diagram

        Returns:
            Problem: a problem with the entities represented explicitly
        """
        self.infer()
        intervals = {}
        for part in self.parts:
            lab = str(self.area2formula(part))
            lab = lab.replace(" ", "").replace("(", "").replace(")", "")
            n = self.parts[part]
            if self.indist[part]:
                e_list = [str(lab)] * n
            else:
                e_list = []
                for j in range(0, n):
                    e = f"{lab}_{j}"
                    e_list.append(e)
            p_int = list2interval(problem, e_list, True)
            intervals[part] = p_int
        for set in self.base_sets:
            area = self.formula2area(set)
            if area in self.parts:
                subsets = [intervals[area]]
            else:
                subsets = [intervals[p] for p in self.get_parts(area)]
            entities = IntervalDict()
            for sub in subsets:
                entities = entities.combine(sub, how=is_distinguishable)
            d = SetFormula(set, entities, problem.universe)
            problem.add_domain(d)
        problem.compute_universe()
        return problem
