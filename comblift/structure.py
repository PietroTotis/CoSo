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

    def __and__(self, rhs):
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
    
    def __str__(self):
        return f"{self.elements}"

    def disjoint(self,rhs):
        inter = self & rhs
        return inter.elements.empty

    def size(self):
        s = 0
        for e in self.elements:
            if not e.empty:
                if e.left == portion.CLOSED and e.right == portion.CLOSED:
                    s += e.upper - e.lower +1
                elif e.left == portion.OPEN and e.right == portion.OPEN:
                    s += e.upper - e.lower -1
                else:
                    s += e.upper - e.lower
        return s

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
        taken = Domain(f"{n}-subset of {self.name}", subset)
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