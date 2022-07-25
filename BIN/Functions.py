#!/usr/bin/env python
"""Summary

Attributes
----------
bin_path : str
    Description
GapSymbol : TYPE
    Description
lib_path : str
    Description
"""
import sys
import os

bin_path = "/lustre/scratch117/cellgen/team297/kt16/BCRSeq/BIN/"
lib_path = "/lustre/scratch117/cellgen/team297/kt16/BCRSeq/LIBRARY/"
if not os.path.exists(bin_path):
    bin_path = os.getcwd() + "/BIN/"
if not os.path.exists(bin_path):
    raise OSError("Cannot locate path to BIN folder")
if not os.path.exists(lib_path):
    lib_path = os.getcwd() + "/LIBRARY/"
if not os.path.exists(lib_path):
    raise OSError("Cannot locate path to LIBRARY folder")
sys.path.append(bin_path)
from collections import defaultdict
import numpy
from numpy import *


class Tree(defaultdict):

    """Summary

    Attributes
    ----------
    value : TYPE
        Description
    """

    def __init__(self, value=None):
        """Summary

        Parameters
        ----------
        value : None, optional
            Description
        """
        super(Tree, self).__init__(Tree)
        self.value = value


def Deconvolute_same_array(tree):
    """Summary

    Parameters
    ----------
    tree : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    decon = Tree()
    inv = {}
    index = 0
    for i in tree:
        if i not in inv:
            index = index + 1
            decon[i][i].value = 1
            inv[i] = i
            array = []
            for j in tree[i]:
                decon[i][j].value = 1
                inv[j] = i
                array.append(j)
            array_new = []
            found = 0
            if len(array) > 0:
                found = 1
            while found == 1:
                array_new = []
                for k in array:
                    for l in tree[k]:
                        if l not in inv:
                            array_new.append(l)
                            inv[l] = i
                            decon[i][l].value = 1
                array = array_new
                if len(array_new) == 0:
                    found = 0
    return (decon, inv)


def Get_codons():
    """Summary

    Returns
    -------
    TYPE
        Description
    """
    fh = open(lib_path + "Codon_table2.txt", "r")
    codon = Tree()
    for l in fh:
        l = l.strip()
        l1 = l.split()
        l2 = list(l1[0])
        codon[l2[0]][l2[1]][l2[2]][l1[1]].value = 1
    fh.close()
    return codon


def Init_rc():
    """Summary

    Returns
    -------
    TYPE
        Description
    """
    b1 = ["A", "T", "G", "C", "N", "."]
    b2 = ["T", "A", "C", "G", "N", "."]
    rc = {}
    for i in range(0, 6):
        rc[b1[i]] = b2[i]
    return rc


def Reverse_comp(seq, rc):
    """Summary

    Parameters
    ----------
    seq : TYPE
        Description
    rc : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    s = ""
    l = len(seq)
    for i in range(0, l):
        j = l - i - 1
        s = s + rc[seq[j]]
    return s


def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp

    Parameters
    ----------
    needle : TYPE
        Description
    haystack : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    m, n = len(needle), len(haystack)
    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m
    row1 = [0] * (n + 1)
    for i in range(0, m):
        row2 = [i + 1]
        for j in range(0, n):
            cost = needle[i] != haystack[j]

            row2.append(
                min(
                    row1[j + 1] + 1,  # deletion
                    row2[j] + 1,  # insertion
                    row1[j] + cost,
                )  # substitution
            )
        row1 = row2
    return min(row1)


import sys
import string


GapSymbol = None


# Utility function.


def hasUniqueElements(L):
    """Summary

    Parameters
    ----------
    L : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    # Does the list have
    # unique elements only?

    d = {}
    for e in L:
        if e == GapSymbol:
            continue
        try:
            f = d[e]
            return 0
        except KeyError:
            d[e] = e
    return 1


class Script:

    """Summary

    Attributes
    ----------
    g : TYPE
        Description
    h : TYPE
        Description
    k : int
        Description
    s : list
        Description
    """

    k = 100
    """Operations transferring one list in another."""

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
        # print
        self.ps("")
        # pass

    def __str__(self):
        """Summary"""
        return str(self.s)

    def __getitem__(self, n):
        """Summary

        Parameters
        ----------
        n : TYPE
            Description
        """
        return self.s[n]

    def __len__(self):
        """Summary"""
        return len(self.s)

    def weight(self, x, y):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
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
        # append "delete k" op
        S = self.s
        if len(S) > 0 and S[-1] != 0:
            # if len(S) > 0 and S[-1] < 0:
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
        # append "insert k" op
        S = self.s

        #   if S[-1] < 0:
        #       S[-1] = k
        #       self.ps("ins.if " + str(k))
        #   else:
        #       S.append(k)
        #       self.ps("ins.else " + str(k))

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
        # append "replace" op
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
        # print script
        # S = self.s
        # print "... %-10s %s" % (note, S)
        pass


class Alignment:

    """Summary"""

    def __init__(self):
        """Summary"""
        # self.GapSymbol = None
        pass

    def diff(self, A, B, M, N, S, tb, te, g, h):
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
                c = c + S.weight(A[1], B[j])
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

        midi = M / 2
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

                c = s + S.weight(A[i], B[j])

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

                c = s + S.weight(A[i + 1], B[j + 1])

                if e < c:
                    c = e
                if d < c:
                    c = d

                s = RR[j]
                RR[j] = c
                SS[j] = d

            S[N] = RR[N]

            # find optimal midpoint

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
                self.diff(A, B, midi, midj, S, tb, g, g, h)
                self.diff(
                    A[midi:], B[midj:], M - midi, N - midj, S, g, te, g, h
                )
            else:
                self.diff(A, B, midi - 1, midj, S, tb, 0.0, g, h)
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
                )

            # return the cost

            return midc

    def do_align(self, A, B, S):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
            Description
        S : TYPE
            Description
        """
        x, y = [], []
        i, j, k, op = 0, 0, 0, 0

        for k in range(len(S)):
            s = S[k]
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

    def align(self, A, B):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
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
        c = self.diff(A, B, M, N, S, G, G, G, H)

        # Removing dummy elements again.
        del A[0]
        del B[0]

        # Create two lists with same length
        # with corresponding matching elements.
        x, y = self.do_align(A, B, S)

        return c, x, y, S

    def partition(self, x, y):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
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


def Do_align(a, b):
    """Summary

    Parameters
    ----------
    a : TYPE
        Description
    b : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    a = list(a)
    b = list(b)
    A = Alignment()
    c, x, y, s = A.align(a, b)
    return (x, y)


def Get_position(a, b):
    """Summary

    Parameters
    ----------
    a : TYPE
        Description
    b : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    a = list(a)
    b = list(b)
    A = Alignment()
    c, x, y, s = A.align(a, b)
    start = -1
    print(y)
    for i in range(0, len(y)):
        if y[i] != None:
            start = i
            break
    return start


"""align.py

Pairwise sequence alignment 
according to Myers & Miller, 1988,
'Optimal alignment in linear space', 
CABIOS, vol. 4, no. 1, pp 11-17.
"""
import os
import sys
import string

GapSymbol = None

# Utility function.


def hasUniqueElements(L):
    """Summary

    Parameters
    ----------
    L : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    # Does the list have
    # unique elements only?

    d = {}
    for e in L:
        if e == GapSymbol:
            continue
        try:
            f = d[e]
            return 0
        except KeyError:
            d[e] = e
    return 1


class Script:

    """Summary

    Attributes
    ----------
    g : TYPE
        Description
    h : TYPE
        Description
    k : int
        Description
    s : list
        Description
    """

    k = 100
    """Operations transferring one list in another."""

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
        # print
        self.ps("")
        # pass

    def __str__(self):
        """Summary"""
        return str(self.s)

    def __getitem__(self, n):
        """Summary

        Parameters
        ----------
        n : TYPE
            Description
        """
        return self.s[n]

    def __len__(self):
        """Summary"""
        return len(self.s)

    def weight(self, x, y):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
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
        # append "delete k" op
        S = self.s
        if len(S) > 0 and S[-1] != 0:
            # if len(S) > 0 and S[-1] < 0:
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
        # append "insert k" op
        S = self.s

        #   if S[-1] < 0:
        #       S[-1] = k
        #       self.ps("ins.if " + str(k))
        #   else:
        #       S.append(k)
        #       self.ps("ins.else " + str(k))

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
        # append "replace" op
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
        # print script
        # S = self.s
        # print "... %-10s %s" % (note, S)
        pass


class Alignment:

    """Summary"""

    def __init__(self):
        """Summary"""
        # self.GapSymbol = None
        pass

    def diff(self, A, B, M, N, S, tb, te, g, h):
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
                c = c + S.weight(A[1], B[j])
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

        midi = M / 2
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

                c = s + S.weight(A[i], B[j])

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

                c = s + S.weight(A[i + 1], B[j + 1])

                if e < c:
                    c = e
                if d < c:
                    c = d

                s = RR[j]
                RR[j] = c
                SS[j] = d

        SS[N] = RR[N]

        # find optimal midpoint

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
            self.diff(A, B, midi, midj, S, tb, g, g, h)
            self.diff(A[midi:], B[midj:], M - midi, N - midj, S, g, te, g, h)
        else:
            self.diff(A, B, midi - 1, midj, S, tb, 0.0, g, h)
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
            )

        # return the cost

        return midc

    def do_align(self, A, B, S):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
            Description
        S : TYPE
            Description
        """
        x, y = [], []
        i, j, k, op = 0, 0, 0, 0

        for k in range(len(S)):
            s = S[k]
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

    def align(self, A, B):
        """Summary

        Parameters
        ----------
        A : TYPE
            Description
        B : TYPE
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
        c = self.diff(A, B, M, N, S, G, G, G, H)

        # Removing dummy elements again.
        del A[0]
        del B[0]

        # Create two lists with same length
        # with corresponding matching elements.
        x, y = self.do_align(A, B, S)

        return c, x, y, S

    def partition(self, x, y):
        """Summary

        Parameters
        ----------
        x : TYPE
            Description
        y : TYPE
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


def printTestCase(A, a, b, x, y, c):
    """Summary

    Parameters
    ----------
    A : TYPE
        Description
    a : TYPE
        Description
    b : TYPE
        Description
    x : TYPE
        Description
    y : TYPE
        Description
    c : TYPE
        Description
    """
    hUE = hasUniqueElements
    print("In-List 1: ", a)
    print("In-List 2: ", b)
    print("Out-List 1:", x)
    print("Out-List 2:", y)
    # print "?", map(hUE, (x, y))
    if list(map(hUE, (x, y))) == [1, 1]:
        print("Partitions:", A.partition(x, y))
    print("Cost:      ", c)


def testCase(note, a, b):
    """Summary

    Parameters
    ----------
    note : TYPE
        Description
    a : TYPE
        Description
    b : TYPE
        Description
    """
    A = Alignment()
    c, x, y, s = A.align(a, b)
    c1, x1, y1, s1 = A.align(b, a)

    print(note)
    print()
    printTestCase(A, a, b, x, y, c)
    if (c, x, y) != (c1, y1, x1):
        print("Asymmetric behaviour:")
        printTestCase(A, b, a, x1, y1, c1)
    print()


def test():
    """Summary

    Raises
    ------
    SystemExit
        Description
    """
    l1 = [1, 2, 3, 4]
    l2 = [1, 2, 1, 4]
    l3 = [1, 2, 3, 4, GapSymbol, GapSymbol]
    res = ["No", "Yes"]

    print("Utility function")
    print()
    for l in [l1, l2, l3]:
        i = hasUniqueElements(l)
        r = res[i]
        print("%s has unique elements: %s" % (l, r))
    print()

    a = []
    b = []
    testCase("Empty lists", a, b)

    a = []
    b = [1, 2, 3]
    testCase("Empty list", a, b)

    a = [1, 2, 3]
    b = [1, 2, 3]
    testCase("Same list", a, b)

    a = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    b = [3, 4, 20, 21, 5, 6, 7, 10, 11]
    testCase("Overlapping", a, b)

    a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    b = [1, 3, 5, 7, 9]
    testCase("Overlapping", a, b)

    a = [1, 3, 5, 7, 9]
    b = [2, 4, 6, 8, 10]
    testCase("Non-overlapping", a, b)

    a = "accgtaccg"
    b = "agtaccccg"
    a = list(a)
    b = list(b)
    testCase("Strings", a, b)

    a = "ac"
    b = "g"
    a = list(a)
    b = list(b)
    testCase("Strings", a, b)

    a = "aaccgtacggt"
    b = "agtgg"
    a = list(a)
    b = list(b)
    testCase("Strings, more complex", a, b)

    # The last case should actually result in:
    #
    #      aaccgtacggt
    #      !---!!--!!-
    #      a   gt  gg
    #
    #      a   gt  gg
    #      !---!!--!!-
    #      aaccgtacggt

    try:
        d1 = sys.argv[1]
        d2 = sys.argv[2]
        a = os.listdir(d1)
        b = os.listdir(d2)
        testCase("Comparing directories %s and %s" % (d1, d2), a, b)
    except IndexError:
        usage = string.split(sys.argv[0], os.sep)[-1]
        msg = "Call %s <dir1> <dir2> " % usage
        msg = msg + "for testing directory compare!"
        print(msg)
        raise SystemExit


def printTestCase2(A, a, b, x, y, c, s):
    """Summary

    Parameters
    ----------
    A : TYPE
        Description
    a : TYPE
        Description
    b : TYPE
        Description
    x : TYPE
        Description
    y : TYPE
        Description
    c : TYPE
        Description
    s : TYPE
        Description
    """
    hUE = hasUniqueElements
    print("In-List 1: ", string.join(a, ""))
    print("In-List 2: ", string.join(b, ""))
    print("Out-List 1:", x)
    print("Out-List 2:", y)
    print("Script:    ", s)
    # print "?", map(hUE, (x, y))
    if list(map(hUE, (x, y))) == [1, 1]:
        print("Partitions:", A.partition(x, y))
    print("Cost:      ", c)


def testCase2(note, a, b):
    """Summary

    Parameters
    ----------
    note : TYPE
        Description
    a : TYPE
        Description
    b : TYPE
        Description
    """
    A = Alignment()
    c, x, y, s = A.align(a, b)
    c1, x1, y1, s1 = A.align(b, a)

    print(note)
    print()
    printTestCase2(A, a, b, x, y, c, s)
    if (c, x, y) != (c1, y1, x1):
        print("Asymmetric behaviour:")
        printTestCase2(A, b, a, x1, y1, c1, s1)
    print()


def test2():
    """Summary"""
    l1 = [1, 2, 3, 4]
    l2 = [1, 2, 1, 4]
    l3 = [1, 2, 3, 4, GapSymbol, GapSymbol]
    res = ["No", "Yes"]

    ##     a = "A very fat cat"
    ##     b = "A very fast cat"
    ##     a = map(None, a)
    ##     b = map(None, b)
    ##     testCase2("Strings", a, b)

    a = "g"
    b = "xx"
    a = list(a)
    b = list(b)
    testCase2("Strings", a, b)


# if __name__ == '__main__':
#  a = "A very fat cat"
#  b = "A very fast cat"
def Do_align(a, b):
    """Summary

    Parameters
    ----------
    a : TYPE
        Description
    b : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    a = list(a)
    b = list(b)
    A = Alignment()
    c, x, y, s = A.align(a, b)
    return (x, y)


def Get_position(a, b):
    """Summary

    Parameters
    ----------
    a : TYPE
        Description
    b : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    a = list(a)
    b = list(b)
    A = Alignment()
    c, x, y, s = A.align(a, b)
    start = -1
    print(y)
    for i in range(0, len(y)):
        if y[i] != None:
            start = i
            break
    return start
