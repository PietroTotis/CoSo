
from util import *


class Venn(object):

    def __init__(self):
        self.base_sets = []
        self.indist = []
        self.sizes = {}
    
    def __repr__(self):
        descr = ""
        for set in self.sizes:
            bsets = [self.base_sets[i] if neg==1 else f"~{self.base_sets[i]}" for i, neg in enumerate(set)]
            descr += f"#({bsets})={self.sizes[set]}\n"
        # descr += str(self.partition)
        return descr

    def add_base_set(self, set, indist):
        self.base_sets.append(set)
        self.indist.append(indist)

    def add_set(self, set, size):
        part = [0] * len(self.base_sets)
        if isinstance(set, str):
            i = self.base_sets.index(set)
            part[i] = 1
            self.sizes[tuple(part)] = size
        else:
            sets = self.flatten_and(set)
            for set in sets:
                if isinstance(set, Not):
                    i = self.base_sets.index(set.child)
                    part[i] = 0
                else:
                    i = self.base_sets.index(set)
                    part[i] = 1
            self.sizes[tuple(part)] = size

    def flatten_and(self, a):
        if isinstance(a,str):
            return [a]
        if isinstance(a.left, str) or isinstance(a.left, Not) :
            return [a.left] + self.flatten_and(a.right)
        elif isinstance(a.right, str) or isinstance(a.right, Not) :
            return [a.right] + self.flatten_and(a.left)
        else:
            return self.flatten_and(a.left) + self.flatten_and(a.right)

        # for s in self.sets:
        #     if self.overlap(set, s):
        #         for c in self.sets[s]:
        #             if self.overlaps(c)

            

         
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
    
#     def get_parts(self, set=None):
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

    
