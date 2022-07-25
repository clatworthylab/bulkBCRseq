#!/usr/bin/env python
"""Summary
"""


class Script:
    """Operations transferring one list in another.

    Attributes
    ----------
    g : TYPE
        Description
    h : TYPE
        Description
    s : list
        Description
    """

    def __init__(self, g=1.5, h=1.0):
        """Summary

        Parameters
        ----------
        g : float, optional
            Description
        h : float, optional
            Description
        """
        self.s = []
        self.g = g
        self.h = h
        self.ps("")

    def __str__(self):
        """Print s.

        Returns
        -------
        TYPE
            Description
        """
        return str(self.s)

    def __getitem__(self, n):
        """Print  s[n].

        Parameters
        ----------
        n : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        return self.s[n]

    def __len__(self):
        """Print len(s).

        Returns
        -------
        TYPE
            Description
        """
        return len(self.s)

    def weight(self, x, y, matrix):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
            Description
        matrix : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        # atomic comparison function
        if x == y:
            return 0
        else:
            return 1

    def gapCost(self, k):
        """Summary

        Parameters
        ----------
        k : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        # k-symbol insert/delete cost
        # return self.g + self.h * k
        if k <= 0:
            return 0
        else:
            return self.g + self.h * k

    def delete(self, k):
        """Summary

        Parameters
        ----------
        k : TYPE
            Description
        """
        S = self.s
        if len(S) > 0 and S[-1] != 0:
            S[-1] = S[-1] - k
            self.ps("del.if " + str(k))
        else:
            S.append(-k)
            self.ps("del.else " + str(k))

    def insert(self, k):
        """Summary

        Parameters
        ----------
        k : TYPE
            Description
        """
        S = self.s
        try:
            if S[-1] < 0:
                S[-1] = k
                self.ps("ins.if " + str(k))
            else:
                S.append(k)
                self.ps("ins.if " + str(k))
        except IndexError:
            S.append(k)
            self.ps("ins.exc " + str(k))

    def replace(self):
        """Summary"""
        S = self.s
        S.append(0)
        self.ps("rep")

    def ps(self, note):
        """Summary

        Parameters
        ----------
        note : TYPE
            Description
        """
        pass


class Alignment:

    """Summary"""

    def __init__(self):
        """Summary"""
        pass

    def diff(self, A, B, M, N, S, tb, te, g, h, matrix):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
            Description
        M : TYPE
            Description
        N : TYPE
            Description
        S : TYPE
            Description
        tb : TYPE
            Description
        te : TYPE
            Description
        g : TYPE
            Description
        h : TYPE
            Description
        matrix : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        # returns the cost of an optimum
        # conversion between A[1..M] and B[1..N]
        # that begins (ends) with a delete
        # if tb (te) is zero and and appends
        # such a conversion to the current script

        # variable setup
        N1 = N + 1
        CC = [0] * N1
        DD = [0] * N1
        RR = [0] * N1
        SS = [0] * N1

        midi, midj, type = 0, 0, 0
        midc = 0.0

        # boundary cases: M <= 1 or N == 0

        if N <= 0 and M > 0:
            S.delete(M)
            return S.gapCost(M)

        if M <= 1:
            if M <= 0:
                S.insert(N)
                return S.gapCost(N)

            if tb > te:
                tb = te

            midc = (tb + h) + S.gapCost(N)
            midj = 0

            for j in range(N + 1):
                c = S.gapCost(j - 1)
                c = c + S.weight(A[1], B[j], matrix)
                c = c + S.gapCost(N - j)
                if c < midc:
                    midc = c
                    midj = j

            if midj == 0:
                S.insert(N)
                S.delete(1)
            else:
                if midj > 1:
                    S.insert(midj - 1)
                S.replace()
                if midj < N:
                    S.insert(N - midj)

            return midc

        # devide: find optimum midpoint
        # (midi, midj) of cost midc

        # Forward phase:
        # Compute C(M/2,k) & D(M/2,k) for all k
        midi = int(M / 2)
        CC[0] = 0.0
        t = g
        for j in range(1, N + 1):
            CC[j] = t = t + h
            DD[j] = t + g

        t = tb
        for i in range(1, midi + 1):
            s = CC[0]
            CC[0] = c = t = t + h
            e = t + g

            for j in range(1, N + 1):
                c = c + g + h
                e = e + h

                if c < e:
                    e = c

                c = CC[j] + g + h
                d = DD[j] + h

                if c < d:
                    d = c

                c = s + S.weight(A[i], B[j], matrix)

                if e < c:
                    c = e
                if d < c:
                    c = d

                s = CC[j]
                CC[j] = c
                DD[j] = d

        DD[0] = CC[0]

        # reverse phase
        # compute R(M/2,k) & S(M/2,k) for all k

        RR[N] = 0.0
        t = g
        for j in range(N - 1, -1, -1):
            RR[j] = t = t + h
            SS[j] = t + g

        t = te
        for i in range(M - 1, midi - 1, -1):
            s = RR[N]
            RR[N] = c = t = t + h
            e = t + g

            for j in range(N - 1, -1, -1):
                c = c + g + h
                e = e + h
                if c < e:
                    e = c

                c = RR[j] + g + h
                d = SS[j] + h

                if c < d:
                    d = c

                c = s + S.weight(A[i + 1], B[j + 1], matrix)

                if e < c:
                    c = e
                if d < c:
                    c = d

                s = RR[j]
                RR[j] = c
                SS[j] = d

        SS[N] = RR[N]

        midc = CC[0] + RR[0]
        midj = 0
        type = 1
        for j in range(0, N + 1):
            c = CC[j] + RR[j]
            if c <= midc:
                if c < midc or CC[j] != DD[j] and RR[j] == SS[j]:
                    midc = c
                    midj = j

        for j in range(N, -1, -1):
            c = DD[j] + SS[j] - g
            if c < midc:
                midc = c
                midj = j
                type = 2

        # conquer: recursively around midpoint
        if type == 1:
            self.diff(A, B, midi, midj, S, tb, g, g, h, matrix)
            self.diff(
                A[midi:], B[midj:], M - midi, N - midj, S, g, te, g, h, matrix
            )
        else:
            self.diff(A, B, midi - 1, midj, S, tb, 0.0, g, h, matrix)
            S.delete(2)
            self.diff(
                A[midi + 1 :],
                B[midj:],
                M - midi - 1,
                N - midj,
                S,
                0.0,
                te,
                g,
                h,
                matrix,
            )

        # return the cost
        return midc

    def do_align(self, A, B, S, GapSymbol=None):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
            Description
        S : TYPE
            Description
        GapSymbol : None, optional
            Description

        Returns
        -------
        TYPE
            Description
        """
        x, y = [], []
        # i, j, k, op = 0, 0, 0, 0
        i, j, k = 0, 0, 0

        for k in range(len(S)):
            s = int(S[k])
            if s == 0:
                if i < len(A) and j < len(B):
                    a, i = A[i], i + 1
                    b, j = B[j], j + 1
                if a == b:
                    x.append(a)
                    y.append(b)
                else:
                    x.append(a)
                    y.append(b)
            elif s < 0:
                y = y + [GapSymbol] * (-s)
                for q in range(i, i - s):
                    x.append(A[q])
                i = i - s
            elif s > 0:
                x = x + [GapSymbol] * (s)
                for q in range(j, j + s):
                    if q < len(B):
                        y.append(B[q])
                    else:
                        y.append(B[q - 1])
                j = j + s

        return x, y

    def align(self, A, B, matrix):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
            Description
        matrix : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        # interface and top level of comparator

        if len(A) == len(B) == 0:
            return 0, A, B, []

        # Build an alignment script.
        G = 1.5
        H = 1.0
        S = Script(G, H)

        # Insert tmp. dummy elements because
        # the algorithm works with indices
        # starting at 1, not 0.
        A.insert(0, 0)
        B.insert(0, 0)
        M = len(A) - 1
        N = len(B) - 1

        # Call recursive function.
        c = self.diff(A, B, M, N, S, G, G, G, H, matrix)

        # Removing dummy elements again.
        del A[0]
        del B[0]

        # Create two lists with same length
        # with corresponding matching elements.
        x, y = self.do_align(A, B, S)

        return c, x, y, S

    def partition(self, x, y, GapSymbol=None):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
            Description
        GapSymbol : None, optional
            Description

        Returns
        -------
        TYPE
            Description
        """
        # assert len(x) == len(y)
        left, both, right = [], [], []

        for i in range(len(x)):
            X, Y = x[i], y[i]
            if X == Y:
                both.append(X)
                continue
            if GapSymbol not in (X, Y):
                left.append(X)
                right.append(Y)
            else:
                if X == GapSymbol:
                    right.append(Y)
                elif Y == GapSymbol:
                    left.append(X)

        return left, both, right
