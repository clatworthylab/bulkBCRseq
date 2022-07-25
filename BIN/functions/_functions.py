#!/usr/bin/env python
"""Summary

Attributes
----------
bin_path : TYPE
    Description
lib_path : TYPE
    Description
main_path : str
    Description
"""
from collections import defaultdict
from ._align import Alignment
import sys
import os

# set up some default paths
main_path = "/lustre/scratch117/cellgen/team297/kt16/BCRSeq/"
bin_path = main_path + "BIN/"
lib_path = main_path + "LIBRARY/"
if not os.path.exists(bin_path):
    bin_path = os.getcwd() + "/BIN/"
    if not os.path.exists(bin_path):
        raise OSError(
            "Cannot locate path to BIN folder. You are currently in {}".format(
                os.getcwd()
            )
        )
    sys.path.append(bin_path)
if not os.path.exists(lib_path):
    lib_path = os.getcwd() + "/LIBRARY/"
    if not os.path.exists(lib_path):
        raise OSError(
            "Cannot locate path to LIBRARY folder. You are currently in {}".format(
                os.getcwd()
            )
        )


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


def fasta_iterator(fh1):
    """Summary

    Parameters
    ----------
    fh1 : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    while True:
        line = fh1.readline()
        if line.startswith(">"):
            break
    while True:
        header = line[1:-1].rstrip()
        sequence = fh1.readline().rstrip()
        while True:
            line = fh1.readline()
            if not line:
                break
            if line.startswith(">"):
                break
            sequence += line.rstrip()
        yield (header, sequence)
        if not line:
            return


def deconvolute_same_array(tree):
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


def get_codons():
    """Summary

    Returns
    -------
    TYPE
        Description
    """
    fh = open(lib_path + "Codon_table2.txt", "r")
    codon = {}
    for l in fh:
        l = l.strip()
        l1 = l.split()
        l2 = list(l1[0])
        codon[l2[0] + l2[1] + l2[2]] = l1[1]
    fh.close()
    return codon


def init_rc():
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


def reverse_comp(seq, rc):
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
        if seq[j] not in rc:
            print(seq)
        else:
            s = s + rc[seq[j]]
    return s


def fuzzy_substring(needle, haystack):
    """
    Calculate the fuzzy match of needle in haystack, using a modified version of the Levenshtein distance algorithm.

    The function is modified from the levenshtein function in the bktree module by Adam Hupp

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
        return needle not in haystack
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


def do_align(a, b):
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


def translate(seq, codon):
    """Summary

    Parameters
    ----------
    seq : TYPE
        Description
    codon : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    p_seq = ""
    for cod in range(0, int(len(seq) / 3 - 1)):
        cod = cod * 3
        c1 = seq[cod + 0] + seq[cod + 1] + seq[cod + 2]
        if c1 in codon:
            p_seq = p_seq + str(codon[c1])
    return p_seq
