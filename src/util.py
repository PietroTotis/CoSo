import portion as P

def is_singleton(interval):
    if interval.left == P.CLOSED and interval.right == P.CLOSED:
        return interval.lower == interval.upper
    if interval.left == P.OPEN and interval.right == P.CLOSED:
        return interval.lower +1 == interval.upper
    if interval.left == P.CLOSED and interval.right == P.OPEN:
        return interval.lower == interval.upper -1
    if interval.left == P.OPEN and interval.right == P.OPEN:
        return interval.lower +1 == interval.upper -1

class Not(object):
    
    def __init__(self, child):
        self.child = child

    def __hash(self):
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

    def __hash(self):
        return hash(str(self))
    
    def __repr__(self):
        l = ("","") if isinstance(self.left, str) else ("(",")")
        r = ("","") if isinstance(self.right, str) else ("(",")")
        return f"{l[0]}{self.left}{l[1]} ∧ {r[0]}{self.right}{r[1]}"

class Or(object):

    def __init__(self, l, r):
        self.left = l
        self.right = r
    
    def __hash(self):
        return hash(str(self))
        
    def __repr__(self):
        l = ("","") if isinstance(self.left, str) else ("(",")")
        r = ("","") if isinstance(self.right, str) else ("(",")")
        return f"{l[0]}{self.left}{l[1]} ∨ {r[0]}{self.right}{r[1]}"