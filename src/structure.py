import portion as P
from util import *


class Domain(object):
    """
    Represents a domain as a portion (set of intervals) with extra info, i.e. distinguishability, size
    and take operation
    
    Attributes
    ----------
    name : str
        the label associated to the domains
    elements : IntervalDict
        intervals corresponding to the entities (mapped to integers) of the domain
        mapped as dict to boolean saying whether elements are distinguishable or not
    """

    def __init__(self, name, elems):
        self.name = name
        self.elements = elems
        self.n_elements = None

    # @staticmethod
    # def is_distinguishable(d1, d2):
    #     # if e is distinguishable truth value in one domain is d1 and the other d2,
    #     # if any of the two elements is distinguishable then keep distinguishing
    #     return d1 or d2 

    # def __and__(self, rhs):
    #     if self.elements in rhs.elements:
    #         i_name = self.name
    #     elif rhs.elements in self.elements:
    #         i_name = rhs.name
    #     else:
    #         i_name = f"({self.name} ∧ {rhs.name})"
    #     i_elem = self.elements & rhs.elements
    #     dist = self.distinguishable.combine(rhs.distinguishable, how=Domain.is_distinguishable)
    #     dist = dist[i_elem]
    #     return Domain(i_name, i_elem, dist)

    def __contains__(self, val):
        return val.elements.domain() in self.elements.domain()

    def __eq__(self, rhs):
        return self.elements == rhs.elements

    def __hash__(self):
        return hash(self.name)

    # def __lt__(self, rhs):
    #     return self in rhs

    # def __or__(self, rhs):
    #     u_name = f"({self.name} ∨ {rhs.name})"
    #     u_elem = self.elements | rhs.elements
    #     dist = self.distinguishable.combine(rhs.distinguishable, how=Domain.is_distinguishable)  
    #     return Domain(u_name, u_elem, dist)
    
    def __sub__(self, rhs):
        c_name = f"({self.name} - {rhs.name})"
        diff =  self.elements - rhs.elements
        dist = self.distinguishable # if distinguishability is different makes no sense
        return Domain(c_name, diff, dist)
        
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        if self.size() > 0 :
            str = f"{self.name} ({self.elements.domain()})"
        else:
            str = f"{self.name} (none)"
        return str

    def all_indistinguishable(self):
        return not (True in self.elements.values())

    def disjoint(self, rhs):
        inter = self & rhs
        return inter.elements.domain().empty

    def size(self):
        if self.n_elements is None:
            s = 0
            for e in self.elements.domain():
                if not e.empty:
                    if e.left == P.CLOSED and e.right == P.CLOSED:
                        s += e.upper - e.lower +1
                    elif e.left == P.OPEN and e.right == P.OPEN:
                        s += e.upper - e.lower -1
                    else:
                        s += e.upper - e.lower
            self.n_elements = s
        return self.n_elements

    def take(self, n):
        """
        returns a subset of itself of size n
        if n > size returns itself
        since we use it for shattering alldiff, assumes that all elements are distinguishable
        """
        iter = P.iterate(self.elements.domain(), step = 1)
        subset = P.empty()
        i = 0
        hasNext = True
        while hasNext and i<n:
            try:
                elem = P.singleton(next(iter))
                subset = subset | elem
            except StopIteration:
                hasNext = False
            i += 1
        dist = P.IntervalDict()
        dist[subset] = True
        taken = Domain(f"{n}x {self.name}", dist)
        return taken

class Structure(object):
    """
    Represents a target structure

    Attributes
    ----------
    name : str
        user-defined name
    type: str
        can be sequence/subset/partition/composition
    spec: bool
        each type has an alternative, if spec is true we use that:
        sequence -> permutation 
        subset -> multisubset
        composition -> multi-composition
        partition -> partition of any size up to n
    domain: DomainFormula
        source set for choices
    size: SizeFormula
        length of sequence/size of subset/number of compositions of partitions    
    """
    def __init__(self, name, type, domain, size = None):
        self.name = name
        self.df = domain
        self.type = str(type)
        self.size = size
        
    def __repr__(self):
        return str(self)

    def __str__(self):
        str = f"{self.type}"
        str += f" ({self.size}) of entity {self.df} ({self.name})" 
        return str

class LiftedSet(object):
    """
    Represents a lifted subset of the universe/a set

    Attributes
    ----------
    name : str
        not important at the moment
    size : SizeFormula
        describes the set of cardinality values that are valid for this set
    cofs : [CountingFormulas]
        describes properties of the set in terms of counting formulas
    relevant : [CountingFormulas]
        defines which groups of indistinguishable objects need to be accounted for
    source : DomainFormula
        the set of which the LiftedSet is subset
    """

    def __init__(self, universe, size, cofs=[]):
        """
        size: portion
        cofs: [CountingFormulas]
        """
        self.name = f"part. of {universe}"
        self.universe = universe
        self.size = size
        self.cofs = self.compact_cofs(cofs)
        self.histogram = {}
        self.check_bound()

    def __and__(self, rhs):
        # name = f"{self.name} /\ {rhs.name}"
        size = self.size & rhs.size
        cofs = self.cofs + rhs.cofs 
        cofs = self.compact_cofs(cofs)
        return LiftedSet(self.universe, size, cofs)
        
    def __contains__(self, constr):
        if hasattr(constr, "formula"):
            for cof in self.cofs:
                if cof.formula == constr.formula:
                    return cof.values in constr.values
            return False
        else:
            return self.size.values in constr.values
        return constr.domain in self.domain

    def __eq__(self, rhs):
        if self.size != rhs.size:
            return False
        elif self.cofs != rhs.cofs:
            return False
        else:
            return True
        
    def __repr__(self):
        return str(self)

    def __str__(self):
        s = f"{self.name}: {self.size}"
        if len(self.cofs) > 0:
            s += "\n"
        s+= "\t counting: "
        s+= ",".join([str(c) for c in self.cofs])
        return s

    def __hash__(self):
        return hash(str(self))
    
    def add_cof(self, cof):
        cofs = self.cofs + [cof]
        return LiftedSet(self.universe, self.size, cofs)

    def bound(self, rv_set):
        bound = False
        if rv_set == self.universe:
            bound = is_singleton(self.size)
        else:
            for cof in self.cofs:
                if cof.formula == rv_set and is_singleton(cof.values):
                    self.histogram[rv_set] = cof.values.lower
                else:
                    bound = bound or is_singleton(cof.values)
        return bound

    def compact_cofs(self, counts):
        compact = []
        remove = []
        updated = False
        indexes = range(0,len(counts))
        for i in indexes:
            cof1 = counts[i]
            for j in [j for j in indexes if j>i]:
                cof2 = counts[j]
                if cof1 != cof2:
                    if cof1.formula == cof2.formula:
                        new_interval = cof1.values & cof2.values
                        merged_cof = cof1.copy()
                        merged_cof.values = new_interval
                        compact.append(merged_cof)
                        remove += [i,j]
                else: # cof1 == cof2:
                    remove += [i]
        keep = [i for i in indexes if not i in remove]
        compact += [counts[i] for i in keep]
        final =  len(remove) == 0
        result = compact if final else self.compact_cofs(compact)
        return result

    def check_bound(self):
        """
        Removes inf upper bound from counting constraints
        """
        for cof in self.cofs:
            if cof.values.upper == P.inf:
                ub = self.size.values.upper +1
                max_int = P.closed(0,ub)
                cof.values = cof.values & max_int

    def copy(self):
        size = self.size.copy()
        cofs = self.cofs.copy()
        return LiftedSet(self.universe, size, cofs)

    def disjoint(self, constr):
        if hasattr(constr, "formula"): # counting formula
            dis_cofs = False
            for cof in self.cofs:
                if cof.formula == constr.formula:
                    if cof.values & constr.values == P.empty():
                        dis_cofs = True
            return dis_cofs
        else: # size formula
            disj = (self.size.values & constr.values) == P.empty()
            return disj

    def feasible(self, rv_set, n):
        if rv_set == self.universe:
            return n in self.size
        else:
            for cof in self.cofs:
                if rv_set in cof.formula:
                    if n>cof.values.upper:
                        return False
                    n_cof = n
                    unfixed = []
                    for rvs in self.histogram:
                        if rvs in cof.formula and rvs!=rv_set:
                            if histogram[rvs]!=-1:
                                unfixed.append(rvs)
                            else:
                                n_cof += histogram[rvs]
                    if len(unfixed) == 0:
                        if n not in cof.values:
                            return False
                    else:
                        sizes = sum([rvs.size() for rvs in unfixed])
                        if sizes + n_cof < cof.values.lower:
                            return False
            fixed = sum([self.histogram[rvs] for rvs in self.histogram if self.histogram[rvs]>-1]) # respect overall size
            if fixed+n not in self.size.values:
                return False
            return True

    def rv_size(self, relevant):
        for cof in self.cofs:
            if cof.formula == relevant:
                return cof.values
                
    def size_is_defined(self):
        s = 0
        for e in self.size.values:
            if not e.empty:
                if e.left == P.CLOSED and e.right == P.CLOSED:
                    s += e.upper - e.lower +1
                elif e.left == P.OPEN and e.right == P.OPEN:
                    s += e.upper - e.lower -1
                else:
                    s += e.upper - e.lower
        return s == 1

    def relevant(self):
        relevant = set([cof.formula for cof in self.cofs])
        if len(relevant) == 0:
            relevant = {self.universe}
        return relevant

    def satisfies(self, constraint):
        sat = None
        for cof in self.cofs:
            if cof.formula in constraint.formula :
                if cof.values in constraint.values:
                    sat = True
                elif constraint.values & cof.values == P.empty:
                    sat = False
        return sat
