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

    Internally areas are identified with tuples of 1/0/-1. Each position in the tuple
    represents a base set, that is, a property. Each position corresponds to a property.
    For each position n, 1 denotes sets of objects that have the n-th property, -1 that 
    they don't have it, 0 that may have it or not.
    """   

    def __init__(self):
        self.assumptions = {}
        self.base_sets = []
        self.base_indist = {}
        self.indist = {}
        self.sets = {}
        self.relevant_areas = []
        self.partitions = {}
        self.universe = None

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
        """
        Add a complex set, i.e. a formula of base sets

        Args:
            set (And/Or/Not/str): a set formula
            size (int): the cardinality of the set
        """
        dnfied = dnfy(set)
        if isinstance(dnfied, Or):
            subsets = [Union(s) for s in flatten(Or, dnfied)]
        else:
            subsets = []
        u = Union(dnfied, size, subsets)
        self.sets[set] = u
        
    def get_subsets(self, a1):
        """
        Check which of the sets with declared size are subsets of set 1

        Args:
            a1 (tuple(int)): the code describing the set

        Returns:
            [tuple(int)]: the codes of the areas contained in a1
        """

        subsets = set()
        for a2 in self.relevant_areas:
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

    def intersects(self, s1, s2):
        """
        Check if two sets are disjoint or share some area

        Args:
            s1 (tuple(int)): area code
            s2 (tuple(int)): area code

        Returns:
            bool: there is an intersection
        """
        for i in range(0, len(s1)):
            if s1[i] == -s2[i] and s1[i] != 0:
                return False
        return True

    def specialize(self, s1, s2):
        """
        Given two sets, generate the codes of the areas that partition s1
        w.r.t. the parts that intersect or not s2

        Args:
            s1 (tuple(int)): set to be specialized
            s2 (tuple(int)): reference set for specialization

        Returns:
            [tuple(int)]: the list of area codes that form s1 w.r.t. the intersections with s2
        """
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

    def get_relevant_areas(self, u):
        """
        Retrieve area codes forming complex sets

        Args:
            u (Union): the set to decompose
        """
        if len(u.subsets) > 0:
            for sub in u.subsets:
                self.get_relevant_areas(sub)
        else:
            code = self.formula2area(u.formula)
            self.relevant_areas.append(code)

    def partition(self, area):
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
        """
        If there is no reference of A/\B in the declaration then assume they are disjoint.
        Update the area code accordingly: add ¬B explicitly to the area code

        Args:
            area (tuple(int)): an area code

        Returns:
            tuple(int): the area code where disjoint subset are explicitly denoted with -1
        """
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
        """
        Builds a Venn area from a complex set. Assumes that we already collected the 
        information about the partitions of the universe composing each set

        Args:
            union (Union): the set

        Returns:
            VennArea: the corresponding VennArea
        """
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
        """
        Create a sort of Venn diagram from the given sets.
        For each set explicitly mentioned we have to identify the atomic parts that compose it, 
        so that later we can do inference on the sizes by summing and subtracting disjoint sets.
        """
        # first collect areas that are mentioned in the program and we have information about their size
        for aset in self.sets.values():
            self.get_relevant_areas(aset)
        # add remaining base sets not explicitly mentioned before
        for aset in self.base_sets:
            self.relevant_areas.append(self.formula2area(aset))
        # then assume disjoint from the information we do not have
        assumed_relevant_areas = []
        for aset in self.relevant_areas:
            ad_area = self.assume_disjoint(aset)
            assumed_relevant_areas.append(ad_area)
            self.assumptions[ad_area] = self.area2formula(aset)
        self.relevant_areas = assumed_relevant_areas
        # Partition each set into disjoint parts
        for area in self.relevant_areas:
            parts = self.partition(area)
            va_parts = [VennArea(self.area2formula(part)) for part in parts]
            self.partitions[self.assumptions[area]] = va_parts
        # Setup universe and add all subsets info
        all_parts = set([part for parts in self.partitions.values() for part in parts])
        self.universe = VennArea("universe", list(all_parts))
        for union in self.sets.values():
            self.universe += self.create_venn_area(union)
        # if there are parts unknown from the set info assume they are empty
        for part in self.universe.parts:
            if part.size is None:
                part.set_size(0)
        # add sets not explicitly mentioned
        for aset in self.base_sets:
            self.universe += self.create_venn_area(Union(aset))
        

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
        self.build_venn()
        intervals = {}
        # create the list of entities from each size
        for part in self.universe.parts:
            if part.indistinguishable:
                e_list = [str(part.name)] * part.size
            else:
                e_list = []
                for j in range(0, part.size):
                    e = f"{part.name}_{j}"
                    e_list.append(e)
            p_int = list2interval(problem, e_list, True)
            intervals[part.name] = p_int
        # create the set formulas for each list of entities
        for subset in self.universe.subsets:
            entities = IntervalDict()
            for part in subset.parts:
                entities = entities.combine(intervals[part.name], how=is_distinguishable)
            d = SetFormula(subset.name, entities, problem.universe)
            problem.add_domain(d)
        problem.compute_universe()
        return problem


class VennArea(object):
    """
    Represents a (sub)set of the universe by keeping track of its atomic parts and subsets
    """

    def __init__(self, name, parts=None, subsets=None, size=None, indistinguishable=True):
        """
        
        Args:
            name (And/Or/Not/str): the set formula corresponding to the set
            parts ([VennArea], optional): The partition of this set. Defaults to [].
            subsets (VennArea, optional): any relevant subset contained in the area. Defaults to [].
            size (int, optional): the cardinality of this set. Defaults to None.
            indistinguishable (bool, optional): whether objects inside are indistinguishable. Defaults to True.
        """
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
        """
        If the two sets have each some exclusive area/objects then combine the two into a parent VennArea
        of which self and other become subsets. Otherwise add the info of the contained set to the parent set

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
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
            return VennArea(f, [va_this_excl, va_inter, va_other_excl], subsets=[self, other])

    def specialize(self, va):
        """
        Add the information of a subset of this area (distinguishability and size)

        Args:
            va (VennArea): an area representing objects contained in self
        """
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
        """
        Split the parts of A and B as A/\B, A/\¬B, ¬A/\B

        Args:
            va (VennArea): other set

        Returns:
            [VennArea]: [A/\B, A/\¬B, ¬A/\B]
        """
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
        """
        Apply the rule that the size of a set is the sum of the sizes of its parts.

        Args:
            n (int): size

        Raises:
            Exception: some counts do not match
        """
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
        """
        Apply the rule that the size of a set is the sum of the sizes of its parts.
        Either one part is missing and we know the other parts' sizes and the partent,
        or we know all parts but not the parent
        """
        for sub in self.subsets:
            sub.infer_size()
        parts_sizes = [p.size for p in self.parts]
        unknown_parts = [p for p in self.parts if p.size is None]
        s = sum([size for size in parts_sizes if size is not None])
        if len(unknown_parts) == 0:
            if self.size is not None and self.size != s:
                raise Exception(f"Inconsistent sizes calculated for {self.name}")
            self.size = s
        if len(unknown_parts) == 1:
            if self.size is not None:
                unknown_parts[0].size = self.size - s
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