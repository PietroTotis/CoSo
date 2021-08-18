
from util import *
from configuration import Domain


# class VennSet(object):
#     # far sparire i Not rimiazzandoli con gli or 

#     def __init__(self, set, size=None):
#         self.formula = set
#         self.size = size
#         self.children = []

#     # negation is always last
#     # assume or formulae are always renamed with children 
#     def contains(self, set):
#         # evaluate formulae
#         if len(self.children) > 0:
#             contains = False
#             for c in self.children:
#                 contains = c.contains(set)        
#                 return contains
#         else:
#             flat = flatten(And, self.formula)
#             if isinstance(set.formula, Or):
#                 # set is in each disjunct
#                 pass
#             else:
#                 flat_set = flatten(And, set.formula)
#                 contains = True
#                 for conjunct in flat:
#                     contains = contains and (conjunct in flat_set)

#     def insert(self, set):
#         # Insert set as And()
#         # Part of the formula can be a child
#         if len(self.children) == 0 and self.contains(set):
#             self.children.append(set)
#         for c in self.children:
#             if c in set:
#                 c.insert(set)

# class Union(object):

#     def __init__(self, formula):
#         self.formula = formula
#         self.susets = []
#         self.size = None

#     # inference: size = sum parts
#     #            #part = size- sum other parts

#     def add_part(self, p):
#         """
#         add a subset to the union, ignoring those already contained in another set
#         """
#         #se due si intersecano cmq considera le tre opzioni
#         # intersezione: somma vettori
#         subset = False
#         for part in self.parts:
#             subset_p = True
#             for i in range(0,len(p)):
#                 ok = p[i] == -1 and part[i] <= 0
#                 ok = p[i] >-1 and p[i]>=part[i]
#                 subset_p = subset_p and ok
#             subset = subset or subset_p
#         if not subset:
#             self.parts.append(p)

class Union(object):

    def __init__(self,formula, size=None, subsets=[]):
        self.formula = formula
        self.size = size
        self.subsets = subsets 

    def __repr__(self):
        descr = f"{self.formula}: {self.size}\n"
        for s in self.subsets:
            descr += f"\t{s}\n"
        return descr

class Venn(object):

    def __init__(self):
        self.base_sets = []
        self.indist = {}
        self.sizes = {}
        self.unions = {}
        self.parts = []
    
    def __repr__(self):
        descr = ""
        for set in self.sizes:
            descr += f"#({self.area2formula(set)}) = {self.sizes[set]}\n"
        # descr += str(self.partition)
        # for u in self.unions:
        #     descr += str(self.unions[u])
        return descr
    
    def add_base_set(self, set, indist):
        self.base_sets.append(set)
        self.indist[set] = indist

    def add_set(self, set, size):
        dnfied = dnfy(set)
        if isinstance(dnfied, Or):
            u = Union(dnfied, size)
            for subset in dnfied:
                u.subsets.append(subset)
            self.unions[set] = u
        else:
            part = self.formula2area(dnfied)
            self.sizes[tuple(part)] = size

    def subsets(self):
        for set in self.sizes:
            subs = self.get_subsets(set)
            if len(subs) > 0:
                name = self.area2formula(set)
                if name not in self.unions:
                    u = Union(name, self.sizes[set], subs)
                    self.unions[name] = u
                else:
                    self.unions[name].subsets += subs
            else:
                self.parts.append(set)
                

    def get_subsets(self, set1):
        """
        Infer from sizes which sets are union of other subsets
        """
        sets = list(self.sizes.keys())
        subsets = set()
        for set2 in sets:
            if set1 != set2:
                if self.is_subset(set2, set1):
                    subsets2 = self.get_subsets(set2)
                    if len(subsets2) > 0:
                        subsets.union(subsets2)
                    else:
                        subsets.add(set2)
                    # print(self.area2formula(set2), "<=", self.area2formula(set1))
        return subsets

    def is_subset(self, set1, set2):
        subset = True
        for i in range(0,len(self.base_sets)):
            if set1[i]!=0 and set1[i] == -set2[i]:
                subset = False
            if set1[i] == 0 and set2[i] == 1:
                subset = False
        return subset

    def infer_union_size(self, u):
        if u.size is None:
            u.size = sum([self.sizes(part) for part in u.subsets])

    def infer_indistinguishability(self):
        for p in self.parts:
            indist = self.sizes[p] > 1 
            sets = [self.base_sets[i] for i, s in enumerate(p) if s==1]
            for set in sets:
                indist = indist and self.indist[set]
            self.indist[p] = indist

    def infer(self):
        self.subsets()
        for u in self.unions:
            print(self.unions[u])
        self.infer_indistinguishability()
        print("parts: ", [self.area2formula(p) for p in self.parts])
        for set in self.indist:
            if isinstance(set, tuple):
                print(self.area2formula(set), self.indist[set])
        for set in self.unions:
            self.infer_union_size(self.unions[set])
        # if u.size != sum(parts) -> phantom rest part
        # u.size = sum(parts)

    def area2formula(self, area):
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
        return part

    # def add_set(self, set, size):
    #     if isinstance(set, str) and not set in self.subsets:
    #         self.subsets[set] = []
    #         self.sizes[set] = size
    #     elif isinstance(set, str) and set in self.subsets:
    #         if self.sizes[set] == size:
    #             pass
    #         else:
    #             print("error")
    #     else:
    #         dnfied = self.dnfy(set)
    #         sets = self.flatten_and(dnfied)
    #         if isinstance(dnfied, Or):
    #             for s in sets:
    #                 if s in self.subsets:
    #                     self.sizes[s] 
    #         else:
                # controlla di quali set è subset
                # trova le foglie
                # aggiungi figlio e (fratelli?)
            

         
# class Venn(object):

#     def __init__(self):
#         self.partition = Partition()
#         self.sets = {}
#         self.sizes = {}

#     def __repr__(self):
#         descr = ""
#         for set in self.sets:
#             descr += f"#{set}={self.sizes[set]} ({self.sets[set]})\n"
#         # descr += str(self.partition)
#         return descr

#     def add_set(self, set, size):
#         base_sets = self.get_base_sets(set)
#         for bs in base_sets:
#             if bs not in self.partition:
#                 self.partition.add_set(bs)
#         parts = self.get_parts(set)
#         self.sets[set] = parts
#         self.sizes[set] = size

#     def get_base_sets(self, set):
#         if isinstance(set, str):
#             return [set]
#         elif isinstance(set, Not):
#             return self.get_base_sets(set.child)
#         else:
#             return self.get_base_sets(set.left)+self.get_base_sets(set.right)

#     def get_parts(self, set):
#         if isinstance(set, Or):
#             return self.get_parts(set.left) + self.get_parts(set.right)
#         elif isinstance(set, And):
#             pleft = self.get_parts(set.left) 
#             pright = self.get_parts(set.right)
#             flat_pleft = []
#             flat_pright = []
#             for part in pleft:
#                 flat_pleft += part.get_parts()
#             for part in pright:
#                 flat_pright += part.get_parts()
#             inter = []
#             for part in flat_pleft:
#                 if part in flat_pright:
#                     inter.append(part)
#             return inter
#         elif isinstance(set, Not):
#             pchild = self.get_parts(set.child)
#             all = self.partition.get_parts()
#             compl =[]
#             for part in all:
#                 if part not in pchild:
#                     compl.append(part)
#             return compl
#         else: # base set
#             return self.partition.get_parts(set) 

# class Partition(object):
#     """
#     Data structure to define a partition tree: level 0 partitions the universe w.r.t. a
#     set S as S and ~S, adding a set adds a new level where the set intersects or not each
#     leaf of the tree.
#     Attributes
#     ----------
#     name : And formula representing the partition as intersection of sets/complements of sets
#     left : name intersect the set corresponding to this level
#     right : name intersect the complement of the set corresponding to this level
#     size : the given size (cardinality) of the set described by name
#     """
#     def __init__(self, part=None, lvl=0):
#         self.part = part
#         self.left = None
#         self.right = None
#         self.size = None
#         self.lvl = lvl

#     def __repr__(self):
#         tabs = "\t"*self.lvl
#         header = tabs + f"{self.part}, size: {self.size}\n"
#         if self .left is not None:
#             return header+f"{str(self.left)} \n {str(self.right)}"
#         else:
#             return header
    
#     def __contains__(self, set):
#         if self.part == set:
#             return True
#         elif self.left is None:
#             return False
#         elif isinstance(self.left.part, And) and self.left.part.left == set:
#             # left starts conjoining set with upper level: it's already in
#             return True
#         else:
#             return set in self.left
    
#     def add_set(self, set):
#         if self.left is None:
#             if self.part is None:
#                 self.left = Partition(set, self.lvl+1)
#                 self.right = Partition(Not(set), self.lvl+1)
#             else:
#                 self.left = Partition(And(set, self.part), self.lvl+1)
#                 self.right = Partition(And(Not(set), self.part), self.lvl+1)
#         else:
#             self.left.add_set(set)
#             self.right.add_set(set)

#     def get_parts(self, set):
#         if isinstance(set, str):
#             return self.get_base_parts(set)
#         elif isinstance(set, Or):
#             lparts = self.get_parts(set.left)
#             rparts = self.get_parts(set.right)
#             lparts = [p for p in lparts if p not in rparts]
#             return lparts+rparts
#         elif isinstance(set, And):
#             lparts = self.get_parts(set.left)
#             rparts = self.get_parts(set.right)
#             return [p for p in lparts if p in rparts]
#         else: # isinstance(set, Not)
#             childparts = self.get_parts(set.child)
#             all = self.get_parts_base()
#             return [p for p in all if p not in childparts]

    
#     def get_parts_base(self, set=None):
#         """
#         set is a base set
#         """
#         if set is None and self.left is None: # if set is None collect all parts
#             return [self.part]
#         elif set is None:
#             return self.left.get_parts() + self.right.get_parts()
#         else:
#             if self.part == set:
#                 return [self]
#             elif self.left is None: # leaf
#                 return []
#             elif isinstance(self.left.part, And) and self.left.part.left == set: # base set
#                 return [self.left]
#             else:
#                 return self.left.get_parts(set) + self.right.get_parts(set)

#     def set_size(self, part, n):
#         if part == self.part:
#             self.size = n
#         else:
#             self.left.set_size(part)
#             self.right.set_size(part)

    
