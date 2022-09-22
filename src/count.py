import math


class Count(object):
    """
    A solution is a set of pairs (count, histogram).
    An histogram is used to keep track of used variables for injectivity constraint.
    For each histogram keep track of the corresponding number of (partial) solutions.
    """

    def __init__(self, count, log, histogram={}, subproblems=1, symbolic=None, tip=""):
        """
        A number plus a bunch of extra info to figure out where it comes from

        Args:
            count (int): the number
            log (ProblemLog): log associated to the count
            histogram (dict, optional): A count for each relevant set used to derive this count. Defaults to {}.
            subproblems (int, optional): Counts the number of subproblems used to derive this count. Defaults to 1.
            symbolic (_type_, optional): A latex representation to display nicely the formula behind this count. Defaults to None.
            tip (str, optional): A text description of the symbolic formula. Defaults to "".
        """
        self.count = count
        self.histograms = [[count, histogram]]
        self.subproblems = subproblems
        self.log = log
        self.symbolic = symbolic if symbolic is not None else str(count)
        self.tip = tip

    def __add__(self, rhs):
        self.count += rhs.count
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
        count = self.count * rhs.count
        subproblems = self.subproblems + rhs.subproblems
        symbolic = self.latex_op("\\cdot", rhs)
        mul = Count(count, self.log, subproblems=subproblems, symbolic=symbolic)
        mul.histograms = histograms
        return mul

    def __sub__(self, rhs):
        if rhs.count == 0:
            return self
        else:
            self.count -= rhs.count
            # self.histograms -= rhs.histograms
            self.subproblems += rhs.subproblems
            self.symbolic = self.latex_op("-", rhs)
            return self

    def __pow__(self, rhs):
        if rhs.count == 1:
            return self
        else:
            self.count = self.count**rhs.count
            self.subproblems += rhs.subproblems
            self.symbolic = self.latex_op("^", rhs)
            return self

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.count} ({self.subproblems} subproblems)"

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
        self.histograms = [[self.count, hst]]

    def with_choices(self, n_choices, tip="Exchangeable choices"):
        updated_hst = []
        for c, hst in self.histograms:
            c = c * n_choices
            updated_hst.append([c, hst])
        self.count *= n_choices
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


class Zero(Count):
    def __init__(self, log, histogram={}):
        super().__init__(0, log, histogram=histogram, subproblems=0, symbolic="")


class One(Count):
    def __init__(self, log, histogram={}):
        super().__init__(1, log, histogram=histogram, subproblems=0, symbolic="")


class Binomial(Count):
    def __init__(self, n, m, log, tip=""):
        count = math.comb(n, m)
        symbolic = f"\\binom{{ {n} }}{{ {m} }}"
        super().__init__(count, log, subproblems=1, symbolic=symbolic, tip=tip)


class Stirling(Count):
    def __init__(self, n, m, log, tip=""):
        count = self.stirling(n, m)
        symbolic = f"\\genfrac{{\\{{}}{{\\}}}}{{0pt}}{{}}{{ {n} }}{{ {m} }}"
        super().__init__(count, log, subproblems=1, symbolic=symbolic, tip=tip)

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
        count = count_n.count // count_m.count
        subproblems = count_n.subproblems + count_m.subproblems
        symbolic = f"\\frac{{{count_n.latex()}}}{{{count_m.latex()}}}"
        super().__init__(
            count,
            count_n.log,
            histogram=histogram,
            subproblems=subproblems,
            symbolic=symbolic,
            tip=tip,
        )


class Factorial(Count):
    def __init__(self, n, log, tip=""):
        count = math.factorial(n)
        symbolic = f"{n}!"
        super().__init__(count, log, subproblems=1, symbolic=symbolic, tip=tip)
