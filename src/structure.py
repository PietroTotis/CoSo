import portion 

def is_singleton(interval):
    if interval.left == portion.CLOSED and interval.right == portion.CLOSED:
        return interval.lower == interval.upper
    if interval.left == portion.OPEN and interval.right == portion.CLOSED:
        return interval.lower +1 == interval.upper
    if interval.left == portion.CLOSED and interval.right == portion.OPEN:
        return interval.lower == interval.upper -1
    if interval.left == portion.OPEN and interval.right == portion.OPEN:
        return interval.lower +1 == interval.upper -1

class Domain(object):
    """
    Represents a domain as a portion (set of intervals)
    
    Attributes
    ----------
    name : str
        the str representation of the formula corresponding to the domains
    elements : Interval
        intervals corresponding to the entities (mapped to integers) of the domain
    distinguishable : IntervalDict
        (sub)sets of (in)distinguishable elements
    """

    def __init__(self, name, elem, distinguishable):
        self.name = name
        self.elements = elem
        self.distinguishable = distinguishable
        self.n_elements = None

    @staticmethod
    def is_distinguishable(d1, d2):
        # if e is distinguishable truth value in one domain is d1 and the other d2,
        # if any of the two elements is distinguishable then keep distinguishing
        return d1 or d2 

    def __and__(self, rhs):
        if self.elements in rhs.elements:
            i_name = self.name
        elif rhs.elements in self.elements:
            i_name = rhs.name
        else:
            i_name = f"({self.name} ∧ {rhs.name})"
        i_elem = self.elements & rhs.elements
        dist = self.distinguishable.combine(rhs.distinguishable, how=Domain.is_distinguishable)
        dist = dist[i_elem]
        return Domain(i_name, i_elem, dist)

    def __contains__(self, val):
        return val.elements in self.elements

    def __eq__(self, rhs):
        return self.elements == rhs.elements

    def __lt__(self, rhs):
        return self in rhs

    def __or__(self, rhs):
        u_name = f"({self.name} ∨ {rhs.name})"
        u_elem = self.elements | rhs.elements
        dist = self.distinguishable.combine(rhs.distinguishable, how=Domain.is_distinguishable)  
        return Domain(u_name, u_elem, dist)
    
    def __sub__(self, rhs):
        c_name = f"¬({rhs.name})"
        diff =  self.elements - rhs.elements
        dist = self.distinguishable[diff]
        return Domain(c_name, diff, dist)
        
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return f"{self.elements}"

    def disjoint(self,rhs):
        inter = self & rhs
        return inter.elements.empty

    def size(self):
        if self.n_elements is None:
            s = 0
            for e in self.elements:
                if not e.empty:
                    if e.left == portion.CLOSED and e.right == portion.CLOSED:
                        s += e.upper - e.lower +1
                    elif e.left == portion.OPEN and e.right == portion.OPEN:
                        s += e.upper - e.lower -1
                    else:
                        s += e.upper - e.lower
            self.n_elements = s
        return self.n_elements

    def take(self, n):
        """
        returns a subset of itself of size n
        if n > size returns itself
        """
        iter = portion.iterate(self.elements, step = 1)
        subset = portion.empty()
        i = 0
        hasNext = True
        while hasNext and i<n:
            try:
                elem = portion.singleton(next(iter))
                subset = subset | elem
            except StopIteration:
                hasNext = False
            i += 1
        taken = Domain(f"{n}x {self.name}", subset)
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
    domain: Constant
        name of the universe, never used actually
    size: Constant
        length of sequence/size of subset/number of compositions of partitions    
    """
    def __init__(self, name, type, spec, domain, size = None):
        self.name = name
        self.domain = domain
        self.type = str(type)
        self.spec = str(spec)=="true"
        self.size = size
        
    def __repr__(self):
        return str(self)

    def __str__(self):
        str = f"{self.size}-{self.type}"
        if self.type == "sequences" and self.spec:
            str = "permutation"
        if self.type == "subset" and self.spec:
            str = f"multi-{str}"
        if self.type == "partition" and self.spec:
            str = f"any-{str}"   
        str += f" of entity {self.domain} ({self.name})" 
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
    source : DomainFormula
        the set of which the LiftedSet is subset
    """

    def __init__(self, name, size, cofs=[]):
        """
        size: portion
        """
        self.name = name
        self.size = size
        self.cofs = self.compact_cofs(cofs)
        self.check_bound()

    def __and__(self, rhs):
        name = f"{self.name} /\ {rhs.name}"
        size = self.size & rhs.size
        cofs = self.cofs + rhs.cofs 
        cofs = self.compact_cofs(cofs)
        return LiftedSet(name, size, cofs)

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
        for c in self.cofs:
            s += f"\t {c}.\n"
        return s

    def __hash__(self):
        return hash(str(self))
    
    def add_cof(self, cof):
        cofs = self.cofs + [cof]
        return LiftedSet(self.name, self.size, cofs)

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
        for cof in self.cofs:
            if cof.values.upper == portion.inf:
                ub = self.size.values.upper +1
                max_int = portion.closed(0,ub)
                cof.values = cof.values & max_int

    def copy(self):
        size = self.size.copy()
        cofs = self.cofs.copy()
        return LiftedSet(self.name, size, cofs)
    
    def size_is_defined(self):
        s = 0
        for e in self.size.values:
            if not e.empty:
                if e.left == portion.CLOSED and e.right == portion.CLOSED:
                    s += e.upper - e.lower +1
                elif e.left == portion.OPEN and e.right == portion.OPEN:
                    s += e.upper - e.lower -1
                else:
                    s += e.upper - e.lower
        return s == 1

    def relevant(self):
        return set([cof.formula for cof in self.cofs])

    def satisfies(self, constraint):
        sat = None
        for cof in self.cofs:
            if cof.formula in constraint.formula :
                if cof.values in constraint.values:
                    sat = True
                elif constraint.values & cof.values == portion.empty:
                    sat = False
        return sat
