import math


class Solution(object):
    def __init__(self, count, log):
        self.count = count
        self.log = log
        self.log.count = count.copy()

    def __str__(self):
        return str(self.count)


class Count(object):
    """
    A solution is a set of pairs (count, histogram).
    An histogram is used to keep track of used variables for injectivity constraint.
    For each histogram keep track of the corresponding number of (partial) solutions.
    """

    def __init__(self, count, histogram={}, subproblems=1, symbolic=None, tip=""):
        """
        A number plus a bunch of extra info to figure out where it comes from

        Args:
            count (int): the number
            histogram (dict, optional): A count for each relevant set used to derive this count. Defaults to {}.
            subproblems (int, optional): Counts the number of subproblems used to derive this count. Defaults to 1.
            symbolic (_type_, optional): A latex representation to display nicely the formula behind this count. Defaults to None.
            tip (str, optional): A text description of the symbolic formula. Defaults to "".
        """
        self.val = count
        self.histograms = [[count, histogram]]
        self.subproblems = subproblems
        self.symbolic = symbolic if symbolic is not None else str(count)
        self.tip = tip

    def __add__(self, rhs):
        self.val += rhs.val
        self.histograms += rhs.histograms
        self.subproblems += rhs.subproblems
        self.symbolic = self.latex_op("+", rhs)
        self.tip = ""
        return self

    def __mul__(self, rhs):
        histograms = []
        for c1, h1 in self.histograms:
            for c2, h2 in rhs.histograms:
                count = c1 * c2
                histogram = self.add_histograms(h1, h2)
                histograms.append([count, histogram])
        count = self.val * rhs.val
        subproblems = self.subproblems + rhs.subproblems
        symbolic = self.latex_op("\\cdot", rhs)
        mul = Count(count, subproblems=subproblems, symbolic=symbolic)
        mul.histograms = histograms
        return mul

    def __sub__(self, rhs):
        if rhs.val == 0:
            return self
        else:
            self.val -= rhs.val
            # self.histograms -= rhs.histograms
            self.subproblems += rhs.subproblems
            self.symbolic = self.latex_op("-", rhs)
            return self

    def __pow__(self, rhs):
        if rhs.val == 1:
            return self
        else:
            self.val = self.val**rhs.val
            self.subproblems += rhs.subproblems
            self.symbolic = self.latex_op("^", rhs)
            return self

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.val} ({self.subproblems} subproblems)"

    def add_histograms(self, h1, h2):
        sum = {}
        for f1 in h1:
            if f1 in h2:
                sum[f1] = h1[f1] + h2[f1]
            else:
                sum[f1] = h1[f1]
        for f2 in h2:
            if f2 not in h1:
                sum[f2] = h2[f2]
        return sum

    def set_histogram(self, hst):
        self.histograms = [[self.val, hst]]

    def with_choices(self, n_choices, tip="Exchangeable choices"):
        updated_hst = []
        for c, hst in self.histograms:
            c = c * n_choices
            updated_hst.append([c, hst])
        self.val *= n_choices
        self.histograms = updated_hst
        self.symbolic = f"{self.add_tip(n_choices, tip)} \cdot {self.symbolic}"
        return self

    def latex(self):
        ignore = self.symbolic == ""
        no_tip = self.tip == ""
        if ignore:
            return ""
        else:
            if no_tip:
                return self.symbolic
            else:
                return self.add_tip(self.symbolic, self.tip)

    def latex_op(self, op, right):
        if self.latex() == "":
            return right.latex()
        elif right.latex() == "":
            return self.latex()
        else:
            return f"{self.latex()} {op} {right.latex()}"

    def add_tip(self, n, tip):
        return f"\\texttip{{ {n} }}{{ {tip} }}"

    def is_zero(self):
        return self.val == 0

    def copy(self):
        return Count(
            self.val, self.histograms, self.subproblems, self.symbolic, self.tip
        )


class Zero(Count):
    def __init__(self, histogram={}, tip="No solution found"):
        super().__init__(0, histogram=histogram, subproblems=0, symbolic="0", tip=tip)

    def __add__(self, rhs):
        return rhs


class One(Count):
    def __init__(self, histogram={}):
        super().__init__(1, histogram=histogram, subproblems=0, symbolic="")


class Binomial(Count):
    def __init__(self, n, m, tip=""):
        count = math.comb(n, m)
        symbolic = f"\\binom{{ {n} }}{{ {m} }}"
        super().__init__(count, subproblems=1, symbolic=symbolic, tip=tip)


class Stirling(Count):
    def __init__(self, n, m, tip=""):
        count = self.stirling(n, m)
        symbolic = f"\\genfrac{{\\{{}}{{\\}}}}{{0pt}}{{}}{{ {n} }}{{ {m} }}"
        super().__init__(count, subproblems=1, symbolic=symbolic, tip=tip)

    def stirling(self, n, k):
        computed = {}

        def stirling_aux(n, k):
            key = str(n) + "," + str(k)
            if key in computed.keys():
                return computed[key]
            if n == k == 0:
                return 1
            if (n > 0 and k == 0) or (n == 0 and k > 0):
                return 0
            if n == k:
                return 1
            if k > n:
                return 0
            result = k * stirling_aux(n - 1, k) + stirling_aux(n - 1, k - 1)
            computed[key] = result
            return result

        return stirling_aux(n, k)


class Divide(Count):
    def __init__(self, count_n, count_m, tip="", histogram={}):
        count = count_n.val // count_m.val
        subproblems = count_n.subproblems + count_m.subproblems
        symbolic = f"\\frac{{{count_n.latex()}}}{{{count_m.latex()}}}"
        super().__init__(
            count,
            histogram=histogram,
            subproblems=subproblems,
            symbolic=symbolic,
            tip=tip,
        )


class Factorial(Count):
    def __init__(self, n, tip=""):
        count = math.factorial(n)
        symbolic = f"{n}!"
        super().__init__(count, subproblems=1, symbolic=symbolic, tip=tip)
