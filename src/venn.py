from portion.dict import IntervalDict
import numpy

from .level_1 import SetFormula
from .util import *


class Union(object):
    """
    Node for a tree-like structure
    """

    def __init__(self, formula, size=None, subsets=[]):
        """
        A node identifies a set with a formula, and registers a size and its subsets

        Args:
            formula (And/Or/Not/str): the formula from base sets
            size (int, optional): the number of objects in the set. Defaults to None.
            subsets (list, optional): the list of subsets. Defaults to [].
        """
        self.formula = formula
        self.size = size
        self.subsets = subsets

    def __repr__(self):
        descr = f"{self.formula}: {self.size}"
        if len(self.subsets) > 0:
            descr += "\n"
            descr += "\n".join([f"\t{s}" for s in self.subsets])
        return descr

    # def set_size(self, s):
    #     if self.atomic():
    #         self.size = s
    #     else:
    #         sum_parts = 0
    #         unknown_parts = []
    #         for p in self.parts:
    #             if p.size is None:
    #                 unknown_parts.append(p)
    #             else:
    #                 sum_parts += p.size
    #         up = len(unknown_parts)
    #         if up == 0:
    #             assert sum_parts == s
    #         elif up == 1:
    #             self.size = s
    #             unknown_parts[0].size = s - sum_parts
    #         else:
    #             self.size = s

    def atomic(self):
        return len(self.subsets) == 0

    # def disjoint(self, aset):
    #     this_inters = aset in self.intersections
    #     part_inters = len(self.part_intersections(aset)) > 0
    #     return this_inters or part_inters

    # def part_intersections(self, aset):
    #     relevant = []
    #     for p in self.parts:
    #         if not p.disjoint(aset):
    #             relevant.append(p)
    #     return relevant

    # def __add__(self, aset):
    #     # 3 cases:
    #     # disjoint - create a new parent
    #     # one in the other - update children
    #     # intersection - new parent, new children
    #     if self.disjoint(aset):
    #         f = Or(self.formula, aset.formula)
    #         if self.size is not None and aset.size is not None:
    #             s = self.size + aset.size
    #         else:
    #             s = None
    #         return Area(f, s, [self, aset])
    #     elif self in aset:
    #         return aset + self
    #     else:
    #         relevant = self.part_intersections(aset)
    #         if len(relevant) == 1:
    #             # one part is either
    #             p = relevant[0]
    #             if p.formula == aset.formula:
    #                 p.size = aset.size
    #                 p.parts = aset.parts
    #             else:

    #         else:


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
        self.relevant_areas = {}
        self.parts = {}
        self.partitions = {}

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
            subsets = [Union(s) for s in flatten(Or, dnfied)]
        else:
            subsets = []
        u = Union(dnfied, size, subsets)
        self.sets[set] = u
        # self.universe += u

    # def get_parts(self, area, f=""):
    #     if f != "" and self.contains(f, area.formula):
    #         return area.parts
    #     else:
    #         parts = []
    #         for p in area.parts:
    #             parts += self.get_parts(p, f)
    #         return parts

    # def explicit_not(self, f):
    #     out = f.child
    #     out_parts = set(self.get_parts(self.universe, out))
    #     all_parts = set(self.get_parts(self.univers))
    #     parts = all_parts - out_parts
    #     formulas = [p.formula for p in parts]
    #     return nest(Or, formulas)

    # def contains(self, f1, f2):
    #     if isinstance(f1, Or):
    #         disjuncts = flatten(Or, f1)
    #         return any([self.contains(d, f2) for d in disjuncts])
    #     elif isinstance(f1, And):
    #         conjuncts = flatten(And, f1)
    #         return all([self.contains(d, f2) for d in conjuncts])
    #     elif isinstance(f1, str):
    #         if isinstance(f2, Or):
    #             disjuncts = flatten(Or, f2)
    #             return all([self.contains(f1, d) for d in disjuncts])
    #         elif isinstance(f2, And):
    #             conjuncts = flatten(And, f2)
    #             return any([self.contains(f1, d) for d in conjuncts])
    #         elif isinstance(f2, str):
    #             return f1 == f2
    #         else:  # Not
    #             f_sets = self.explicit_not(f2)
    #             return self.contains(f1, f_sets)
    #     else:
    #         f_sets = self.explicit_not(f1)
    #         return self.contains(f2, f_sets)

    # def subsets(self):
    #     """
    #     For all base sets collect subsets from declared intersections
    #     """
    #     for set in self.base_sets:
    #         subs = self.get_subsets(set)
    #         if set in self.sets and len(subs) > 0:
    #             s = self.sets[set].size
    #             self.sets[set] = Union(set, s, subs)
    #         elif set in self.sets and len(subs) > 0:
    #             s = None
    #             self.sets[set].subsets += subs
    #         elif set not in self.sets and len(subs) == 0:
    #             self.sets[set] = Union(set)
    #         else:
    #             pass
    #             # self.parts.append(set)

    def get_subsets(self, a1):
        """
        Check which of the sets with declared size are subsets of set 1
        """

        subsets = set()
        # a1 = self.formula2area(set1)
        all_areas = [a for areas in self.relevant_areas.values() for a in areas]
        for a2 in all_areas:
            # a2 = self.formula2area(set2)
            if a1 != a2:
                if self.is_subset(a2, a1):
                    subsets2 = self.get_subsets(a2)
                    if len(subsets2) > 0:
                        subsets.union(subsets2)
                    else:
                        subsets.add(a2)
        return list(subsets)

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
            if s1[i] == 0 and s2[i] != 0:
                subset = False
        return subset

    # def add_parts(self, set1, set2):
    #     a1 = self.formula2area(set1.name)
    #     a2 = self.formula2area(set2.name)
    #     s1 = set1.size
    #     s2 = set2.size
    #     for sub in set1:
    #         if set2 in sub:
    #             self.add_parts(sub, set2)
    #     if s1 is not None and s2 is not None:

    #         if s1 == s2:

    #     # size defined, equal
    #     # size defined, different
    #     if set1.size == set2.size:

    def intersects(self, s1, s2):
        for i in range(0, len(s1)):
            if s1[i] == -s2[i] and s1[i] != 0:
                return False
        return True

    def specialize(self, s1, s2):
        parts = [[]]
        for i, c1 in enumerate(s1):
            c2 = s2[i]
            if c1 == c2:
                for p in parts:
                    p.append(c1)
            else:
                part1 = [p + [1] for p in parts]
                part2 = [p + [-1] for p in parts]
                parts = part1 + part2
        return parts

    def add_relevant_area(self, subset):
        code = self.formula2area(subset)
        level = len([i for i in code if i != 0])
        if level not in self.relevant_areas:
            self.relevant_areas[level] = {code}
        else:
            self.relevant_areas[level] |= {code}

    def setup_area(self, u):
        if len(u.subsets) > 0:
            for sub in u.subsets:
                self.setup_area(sub)
        else:
            self.add_relevant_area(u.formula)

    def setup_areas(self):
        for u in self.sets.values():
            self.setup_area(u)

    def partition(self, area):
        # aset = self.area2formula(area)
        # if isinstance(aset, Or):
        #     subparts = {
        #         part for sub in flatten(Or, aset) for part in self.partition(sub)
        #     }
        #     return subparts
        # else:
        parts = []
        subs = self.get_subsets(area)
        for i, s1 in enumerate(subs):
            intersections = [s2 for s2 in subs[i + 1 :] if self.intersects(s1, s2)]
            if len(intersections) == 0:
                parts.append(s1)
            else:
                for s2 in intersections:
                    parts += self.specialize(s1, s2)
        rest = list(area)
        for i in range(0,len(area)):
            if area[i] == 0 and any([p[i]==1 for p in parts]):
                rest[i] = -1        
        parts.append(tuple(rest))
        return parts

    def assume_disjoint(self, area):
        subs = self.get_subsets(area)
        # count declared intersections with any base set
        m = tuple(map(lambda x: x.count(1), zip(*subs)))
        if -1 in area:
            # ignore complements: disjoint assumption not applicable
            return area
        assumed = []
        if len(m) == 0:
            # no subsets: already a part, assume disjoint all non-explicit (1) base sets
            assumed = [a if a == 1 else -1 for a in area]
        else:
            for i, a in enumerate(area):
                if a != 0:
                    # declared
                    assumed.append(a)
                else:
                    if m[i] > 0:
                        # there are intersections: both -1 and 1 are possible
                        assumed.append(0)
                    else:
                        # no other intersections: assume disjoint
                        assumed.append(-1)
        return tuple(assumed)

    def create_venn_area(self, union):
        if len(union.subsets) >0:
            subs = []
            for sub in union.subsets:
                subs.append(self.create_venn_area(sub))
            va = subs[0]
            for s in subs[1:]:
                va += s
            va.set_size(union.size)
            return va
        else:
            parts = self.partitions[union.formula]
            indist = self.base_indist.get(union.formula, True)
            return VennArea(union.formula, parts, size=union.size, indistinguishable=indist)


    def build_venn(self):
        self.setup_areas()
        for i in self.relevant_areas:
            self.relevant_areas[i] = [self.assume_disjoint(area) for area in self.relevant_areas[i]]
        # build all non-empty partitions
        for i in self.relevant_areas:
            for aset in self.relevant_areas[i]:
                parts = [VennArea(self.area2formula(part)) for part in self.partition(aset)]
                self.partitions[self.area2formula(aset)] = parts
        #         print(i, parts)
        all_parts = set([part for parts in self.partitions.values() for part in parts])
        # univ_parts = [VennArea(part, []) for part in all_parts]
        # univ = (0,)*len(self.base_sets)
        universe = VennArea("universe", list(all_parts))
        print(universe)
        for union in self.sets.values():
            universe += self.create_venn_area(union)
        print(universe)
        print("fine")
        # for aset in self.partitions:
        #     print(aset)
        #     parts = [VennArea(p, []) for p in self.partitions[aset]]
        #     va = VennArea(aset, parts)
        #     print(f"adding {aset}, {parts}") 
        #     universe += va
        # print(universe)


        
        # for si in self.partitions:
        #     print(
        #         self.area2formula(si),
        #         [self.area2formula(s) for s in self.partitions[si]],
        #     )
        # add sets as union of partitions

        # univ_area = (0,) * len(self.base_sets)
        # universe = VennArea(univ_area, nonempty_sets, name="universe")
        # for s in self.base_sets:
        #     vas = VennArea(self.formula2area(s), s)
        #     universe.add_set(vas)
        #     print(universe)
        # parts = self.get_parts(univ_area)
        # for p in parts:
        #     f = self.area2formula(p)
        #     if f in self.sets:
        #         self.parts[p] = self.sets[f].size
        # return parts

    # def infer_size(self, parts):
    #     if isinstance(subset, Union):
    #         self.infer_union_size(subset)
    #     else:  # tuple
    #         if subset in self.parts:
    #             return self.parts[subset]
    #         elif 0 not in subset:
    #             return 0
    #         else:
    #             return 0

    # def infer_union_size(self, u):
    #     """
    #     Inference on the sizes of the sets

    #     Args:
    #         u (Union): a set

    #     Raises:
    #         Exception: we might find inconsistent sizes
    #         Exception: we might miss some sizes
    #     """
    #     f = self.formula2area(u.formula)
    #     parts = [tuple(p) for p in self.get_parts(f)]
    #     known = [p for p in parts if p in self.parts]
    #     ksize = sum([self.parts[part] for part in known])
    #     unknown = [p for p in parts if p not in self.parts]
    #     if len(unknown) == 0:
    #         # the set size is the sum of the parts' sizes
    #         if u.size is None:
    #             u.size = ksize
    #         elif u.size != ksize:
    #             raise Exception(
    #                 f"Inconsistent size of {u.formula} (expected {u.size}, got {ksize})"
    #             )
    #     elif len(unknown) == 1:
    #         # if the size is not given, assume it is sum of known parts
    #         # otherwise the size of unknown is the difference
    #         if u.size is None:
    #             u.size = ksize
    #             self.parts[unknown[0]] = 0
    #         else:
    #             self.parts[unknown[0]] = u.size - ksize
    #     else:
    #         # many unknown: try to infer from subsets assuming they are disjoint if
    #         # not given otherwise
    #         u.size = 0
    #         for i, s1 in enumerate(u.subsets):
    #             # print(s1, self.area2formula(s1))
    #             u.size += self.infer_size(s1)
    #             # for s2 in u.subsets[i+1:]:
    #             #     u.size -= self.infer_union_size(s1&s2)
    #     return u.size
    # if u.size is None:
    #     if len(known) == 0:
    #         raise Exception(
    #             f"I don't know neither the size of {u.formula} nor any of its subsets"
    #         )
    #     else:
    #         u.size = ksize
    # else:
    #     if len(known) == 0:
    #         self.parts[f] = u.size
    # missing the case where both #A and #B are unknown but #A u B is known

    def infer_indistinguishability(self):
        for p in self.parts:
            indist = self.parts[p] > 1
            sets = [self.base_sets[i] for i, s in enumerate(p) if s == 1]
            for set in sets:
                indist = indist and self.base_indist[set]
            self.indist[p] = indist

    def infer(self):
        # subset_tree = self.subsets()
        self.build_venn()
        # for set in self.sets:
        #     self.infer_size(self.sets[set])
        # self.infer_indistinguishability()

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


class VennArea(object):

    def __init__(self, name, parts=None, subsets=None, size=None, indistinguishable=True):
        self.name = name
        self.parts = [] if parts is None else parts
        self.subsets = [] if subsets is None else subsets
        self.size = size
        self.indistinguishable = indistinguishable

    def __hash__(self):
        if len(self.parts)>0:
            return hash(tuple(self.parts))
        else:
            return hash(self.name)

    def __eq__(self, other):
        # print(self.name, self.parts)
        if len(self.parts) == 0:
            return self.name == other.name
        else:
            return self in other and other in self

    def __add__(self, other):
        if self == other:
            self.set_size(other.size)
            return self    
        elif self in other:
            other.specialize(self)
            return other
        elif other in self:
            self.specialize(other)
            return self
        else:   
            dist = not(self.indistinguishable and other.indistinguishable)
            this_excl, inter, other_excl = self.intersect(other)
            f_this_excl = nest(Or, [p.name for p in this_excl])
            va_this_excl = VennArea(f_this_excl, parts=this_excl, indistinguishable=self.indistinguishable)
            f_inter = nest(Or, [p.name for p in inter])
            va_inter = VennArea(f_inter, parts=inter, indistinguishable=dist)
            f_other_excl = nest(Or, [p.name for p in other_excl])
            va_other_excl = VennArea(f_other_excl, parts=other_excl, indistinguishable=other.indistinguishable)
            f = nest(Or,[f_this_excl, f_inter, f_other_excl])
            return VennArea(f, [va_this_excl, va_inter, va_other_excl])

    def specialize(self, va):
        _, inter, _ = self.intersect(va)
        if not va.indistinguishable:
            for i in inter:
                i.set_distinguishable()
        if len(inter) >1:
            va.parts = inter
            self.subsets.append(va)
        else:
            if inter[0].size is None:
                inter[0].size = va.size
        self.infer_size()

    def intersect(self, va):
        this_exclusive = [p for p in self.parts if p not in va.parts]
        intersection = [p for p in self.parts if p in va.parts]
        va_exclusive = [p for p in va.parts if p not in self.parts]
        return this_exclusive, intersection, va_exclusive
    
    def intersects(self, va):
        _, inter, _ = self.intersect(va)
        return len(inter)

    def __contains__(self, va):
        _, _, va_excl = self.intersect(va)
        return len(va_excl) == 0

    def name(self):
        return self.names[0]
    
    def set_size(self, n):
        if n is None:
            return
        if self.size is None:
            self.size = n
            part_sizes = [p.size for p in self.parts]
            unknown_parts = [p for p in self.parts if p.size is None]
            if len(unknown_parts)==1:
                sum_known = sum([p for p in part_sizes if p is not None])
                unknown_parts[0].set_size(n-sum_known)
        else:
            self.infer_size()
            if self.size != n:
                raise Exception(f"Inconsistent sizes for {self.name}")
            
    def set_distinguishable(self):
        self.indistinguishable = False
        for sub in self.subsets:
            sub.set_distinguishable()
        for part in self.parts:
            part.set_distinguishable()

    def infer_size(self):
        for sub in self.subsets:
            sub.infer_size()
        parts_sizes = [p.size for p in self.parts]
        if self.name == "red":
            print(self.parts, parts_sizes)
        unknown_parts = [p for p in self.parts if p.size is None]
        s = sum([size for size in parts_sizes if size is not None])
        if len(unknown_parts) == 0:
            if self.size is not None and self.size != s:
                raise Exception(f"Inconsistent sizes calculated for {self.name}")
            self.size = s
            print(f"inferred {self.name}")
        if len(unknown_parts) == 1:
            if self.size is not None:
                unknown_parts[0].size = self.size - s
                print(f"inferred {unknown_parts[0].name}")
        if self.size == 1:
            self.indistinguishable = False


    def pretty_print(self, indent):
        space = indent*"  "
        name = str(self.name) if not self.indistinguishable else f"indist. {self.name}"
        s = f"{space}{name} ({self.size} obj)"
        s += f" > " + ", ".join([f"{sub.name} ({sub.size})" for sub in self.subsets])
        for p in self.parts:
            s += f"\n{space}{p.pretty_print(indent+1)}"
        return s


    def __repr__(self):
        return self.pretty_print(0)