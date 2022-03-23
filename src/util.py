import portion as P
import os
from portion.dict import IntervalDict

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))

def is_distinguishable(d1, d2):
    # if e is distinguishable truth value in one domain is d1 and the other d2,
    # if any of the two elements is distinguishable then keep distinguishing
    return d1 or d2 

def interval_closed(interval, lb_default=0, ub_default=P.inf):
    if interval.lower != P.inf and interval.upper!=P.inf:
        lb = interval.lower if interval.left == P.CLOSED else interval.lower +1
        ub = interval.upper if interval.right == P.CLOSED else interval.upper -1
        return (lb,ub)
    elif interval.lower != P.inf:
        lb = interval.lower if interval.left == P.CLOSED else interval.lower +1
        return (lb,ub_default)
    else:
        ub = interval.upper if interval.right == P.CLOSED else interval.upper -1
        return (lb_default,ub)

def is_singleton(interval):
    if interval.lower == P.inf or interval.upper==P.inf:
        return False
    if interval.left == P.CLOSED and interval.right == P.CLOSED:
        return interval.lower == interval.upper
    if interval.left == P.OPEN and interval.right == P.CLOSED:
        return interval.lower +1 == interval.upper
    if interval.left == P.CLOSED and interval.right == P.OPEN:
        return interval.lower == interval.upper -1
    if interval.left == P.OPEN and interval.right == P.OPEN:
        return interval.lower +1 == interval.upper -1

def combine(l, r):
    # update indistinguishability when and/or domains
    s_dom = l.elements.domain()
    r_int = r.elements.domain()
    if not s_dom.overlaps(r_int):
        return l + r
    else:
        comb = IntervalDict()
        l_keys = l.elements.keys()
        r_keys = r.elements.keys()
        # print(l,"//", r, "   ")
        l_n = 0
        r_n = 0
        while l_n < len(l_keys) or r_n < len(r_keys):
            # iterate partitions on left/right 
            # print(l_n, "----", r_n )
            if l_n == len(l_keys):  # add non-overlapping left 
                r_int = r_keys[r_n]
                comb = comb.combine(r.elements[r_int], how=is_distinguishable)
                r_n += 1
            elif r_n == len(r_keys): # add non-overlapping right 
                l_int = l_keys[l_n]
                comb = comb.combine(l.elements[l_int], how=is_distinguishable)
                l_n += 1
            else: # left and right partition share elements
                l_int = l_keys[l_n]
                r_int = r_keys[r_n]
                # print("l:", l_int, "r:", r_int)
                if l_int<r_int:
                    comb = comb.combine(l.elements[l_int], how=is_distinguishable)
                    l_n += 1
                elif r_int<l_int:
                    comb = comb.combine(r.elements[r_int], how=is_distinguishable)
                    comb[r_int] = r.elements[r_int]
                    r_n += 1
                else:
                    # print(l_int,"--",r_int)
                    l_out = l_int - r_int
                    inter = l_int & r_int
                    r_out = r_int - l_int
                    # print("vals", l_out,inter,r_out)
                    if not l_out.empty:
                        if is_singleton(l_out):
                            comb[l_out] = True
                        else: # partiton of size 1 is always distinguishable
                            comb = comb.combine(l.elements[l_out], how=is_distinguishable)
                    if not r_out.empty:
                        if is_singleton(r_out):
                            comb[r_out] = True
                        else:
                            comb = comb.combine(r.elements[r_out], how=is_distinguishable)
                    if not inter.empty:
                        if is_singleton(inter):
                            comb[inter] = True
                        else:
                            comb = comb.combine(l.elements[inter], how=is_distinguishable)
                            comb = comb.combine(r.elements[inter], how=is_distinguishable)
                    l_n += 1
                    r_n += 1
        return comb


class Not(object):
    
    def __init__(self, child):
        self.child = child

    def __eq__(self, rhs):
        if isinstance(rhs, Not):
            return self.child==rhs.child
        else:
            return False

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        if isinstance(self.child,str):
            return f"¬{self.child}"
        else:
            return f"¬({self.child})"

class And(object):

    def __init__(self, l, r):
        self.left = l
        self.right = r

    def __eq__(self, rhs):
        if isinstance(rhs, And):
            same = self.left == rhs.left and self.right == rhs.right
            inverted = self.right == rhs.left and self.left == rhs.right
            return same or inverted
        else:
            return False

    def __hash__(self):
        return hash(str(self))
    
    def __repr__(self):
        # if self.left.name == "":
        #     return str(self.right)
        # elif self.right.name == "":
        #     return str(self.left)
        # else:
        l = ("","") if isinstance(self.left, str) else ("(",")")
        r = ("","") if isinstance(self.right, str) else ("(",")")
        return f"{l[0]}{self.left}{l[1]} ∧ {r[0]}{self.right}{r[1]}"
    
    

class Or(object):

    def __init__(self, l, r):
        self.left = l
        self.right = r
    
    def __hash__(self):
        return hash(str(self))
        
    def __eq__(self, rhs):
        if isinstance(rhs, Or):
            same = self.left == rhs.left and self.right == rhs.right
            inverted = self.right == rhs.left and self.left == rhs.right
            return same or inverted
        else:
            return False

    def __repr__(self):
        l = ("","") if isinstance(self.left, str) else ("(",")")
        r = ("","") if isinstance(self.right, str) else ("(",")")
        return f"{l[0]}{self.left}{l[1]} ∨ {r[0]}{self.right}{r[1]}"

def dnfy(set):
        """
        Transform a set formula into dnf
        """
        if isinstance(set, str):
            return set
        if isinstance(set, Or):
            ldnfy = dnfy(set.left)
            rdnfy = dnfy(set.right)
            return Or(ldnfy, rdnfy)
        elif isinstance(set, Not):
            if isinstance(set.child, Not):
                return dnfy(set.child.child)
            elif isinstance(set.child, And):
                ldnfy = dnfy(Not(set.child.left))
                rdnfy = dnfy(Not(set.child.right))
                return Or(ldnfy,rdnfy)
            elif isinstance(set.child, Or):
                ldnfy = dnfy(Not(set.child.left))
                rdnfy = dnfy(Not(set.child.right))
                return And(ldnfy,rdnfy)
            else:
                return set
        else: #isinstance(set, And):
            if isinstance(set.right, Or):
                rdnfy = dnfy(set.right)
                return Or(And(set.right.left, rdnfy), And(set.right.left, rdnfy))
            elif isinstance(set.right,Or):
                ldnfy = dnfy(set.left)
                return Or(And(set.right.left, ldnfy), And(set.right.left, ldnfy))
            else:
                return set

def flatten(op, a):
    """
    Flatten conjunction or disjunction from binary tree to list
    op : And/Or
    a : And/Or/Not/str
    """
    if isinstance(a, str) or isinstance(a, Not):
        return [a]
    if isinstance(a, op) and (isinstance(a.left, str) or isinstance(a.left, Not)):
        return [a.left] + flatten(op, a.right)
    elif isinstance(a, op) and (isinstance(a.right, str) or isinstance(a.right, Not)):
        return [a.right] + flatten(op, a.left)
    elif isinstance(a, op) and isinstance(a.left, op) and isinstance(a.right, op):
        return flatten(op, a.left) + flatten(op, a.right)
    else: # flattening And (Or) but both children are Or (And)
        return a 

def nest(op, l):
    """
    Turn a list l into a conjunction or disjunction binary tree
    op : And/Or
    l : list of str
    """
    if len(l)==1:
        return l[0]
    else:
        return op(l[0], nest(op, l[1:]))