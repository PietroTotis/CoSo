import portion 


class Domain(object):
    """
    Represents a domain as a portion (set of intervals)
    
    Attributes
    ----------
    name : str
        the str representation of the formula corresponding to the domains
    elements : portion
        intervals corresponding to the entities (mapped to integers) of the domain
    """

    def __init__(self, name, elem):
        self.name = name
        self.elements = elem
        self.n_elements = None

    def __and__(self, rhs):
        if self.elements in rhs.elements:
            i_name = self.name
        elif rhs.elements in self.elements:
            i_name = rhs.name
        else:
            i_name = f"({self.name} ∧ {rhs.name})"
        i_elem = self.elements & rhs.elements
        return Domain(i_name, i_elem)

    def __contains__(self, val):
        return val.elements in self.elements

    def __eq__(self, rhs):
        return self.elements == rhs.elements

    def __lt__(self, rhs):
        return self in rhs

    def __or__(self, rhs):
        u_name = f"({self.name} ∨ {rhs.name})"
        u_elem = self.elements | rhs.elements
        return Domain(u_name, u_elem)
    
    def __sub__(self, rhs):
        c_name = f"¬({rhs.name})"
        return Domain(c_name, self.elements - rhs.elements)
        
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
    constraints : [CountingFormulas]
        describes properties of the set in terms of counting formulas
    source : DomainFormula
        the set of which the LiftedSet is subset
    """

    def __init__(self, name, size, source, constraints=[]):
        """
        size: portion
        """
        self.name = name
        self.size = size
        self.source = source
        self.constraints = constraints

    def __and__(self, rhs):
        name = f"{self.name} /\ {rhs.name}"
        size = self.size & rhs.size
        constraints = self.constraints + rhs.constraints #compact
        return LiftedSet(name, size, self.source, constraints)

    def __eq__(self, rhs):
        if self.size != rhs.size:
            return False
        elif self.constraints != rhs.constraints:
            return False
        else:
            return True
        
    def __repr__(self):
        return str(self)

    def __str__(self):
        s = f"{self.source}: {self.size}"
        if len(self.constraints) > 0:
            s += "\n"
        for c in self.constraints:
            s += f"\t {c}.\n"
        return s

    def __hash__(self):
        return hash(str(self))
    
    def copy(self):
        new = LiftedSet(self.name, self.size, self.source)
        new.constraints = self.constraints
        return new
    
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
        return set([cof.formula for cof in self.constraints])

    def satisfies(self, constraint):
        sat = None
        for cof in self.constraints:
            if cof.formula in constraint.formula :
                if constraint.values in cof.values:
                    sat = True
                elif constraint.values & cof.values == portion.empty:
                    sat = False
        return sat

