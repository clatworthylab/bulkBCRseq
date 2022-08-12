#!/usr/bin/env python
import re
import subprocess

from collections import defaultdict
from pathlib import Path
from typing import Tuple, Dict, List

from isotyper.utilities._settings import (
    REVERSE_COMPLEMENT_DICT,
    READ_NUMBER_DIVISION,
)


class Tree(defaultdict):
    """Tree data structure for nested dictionaries."""

    def __init__(self, value=None):
        """Initialise tree.

        Parameters
        ----------
        value : None, optional
            If specified, creates an entry e.g.
            Tree()['key']['entry'].value = 1
        """
        super(Tree, self).__init__(Tree)
        self.value = value


def fasta_iterator(fh1: Path) -> Tuple[str, str]:
    """Fasta reader.

    Parameters
    ----------
    fh1 : Path
        path to fasta file.
    Returns
    -------
    Tuple[str, str]
        pair of header and sequence from a fasta file.
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


def deconvolute_same_array(tree: Tree) -> Tuple[Tree, Dict]:
    """Summary
    Parameters
    ----------
    tree : Tree
        Tree data structure
    Returns
    -------
    Tuple[Tree, Dict]
        Tree data and dictionary.
    """
    decon = Tree()
    inv = {}
    index = 0
    for i in tree:
        if i not in inv:
            index += 1
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


def get_codons(codon_file: Path) -> Dict:
    """Get codons from codon table file.

    Parameters
    ----------
    codon_file : Path
        path to codon file e.g. Codon_table2.txt

    Returns
    -------
    Dict
        Codons in a dictionary.
    """
    fh = open(codon_file, "r")
    codon = {}
    for l in fh:
        l = l.strip()
        l1 = l.split()
        l2 = list(l1[0])
        codon[l2[0] + l2[1] + l2[2]] = l1[1]
    fh.close()
    return codon


def reverse_comp(seq: str) -> str:
    """convert nucleotide sequence to reverse complementary.
    Parameters
    ----------
    seq : str
        nucleotide sequence
    Returns
    -------
    str
        reverse complementary sequence
    """
    s = ""
    l = len(seq)
    for i in range(0, l):
        j = l - i - 1
        if seq[j] in REVERSE_COMPLEMENT_DICT:
            s = s + REVERSE_COMPLEMENT_DICT[seq[j]]
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


def translate(seq: str, codon: Dict) -> str:
    """Translate sequence to amino acid.

    Parameters
    ----------
    seq : str
        sequence to translate.
    codon : Dict
        dictionary of codons.
    Returns
    -------
    str
        translated sequence.
    """
    p_seq = ""
    for cod in range(0, int(len(seq) / 3 - 1)):
        cod = cod * 3
        c1 = seq[cod + 0] + seq[cod + 1] + seq[cod + 2]
        if c1 in codon:
            p_seq = p_seq + str(codon[c1])
    return p_seq


def write_out(out: str, file: Path):
    """Generic writer in append mode.

    Parameters
    ----------
    out : str
        string to write.
    file : Path
        path of file to write.
    """
    fh = open(file, "a")
    fh.write(out)
    fh.close()


def create_file(file: Path):
    """Create a file.

    Parameters
    ----------
    file : Path
        path of file to create.
    """
    fh = open(file, "w")
    fh.close()


def join_reads(s1: str, s2: str, length: int) -> Tuple[str, int]:
    """Join reads.

    Parameters
    ----------
    s1 : str
        sequence 1
    s2 : str
        sequence 2
    length : int
        length to trim

    Returns
    -------
    Tuple[str, int]
        joint read and log of whether it failed.
    """
    seq = ""
    (s2) = reverse_comp(s2)
    (l1, l2) = (len(s1), len(s2))
    failed = 1
    for i in range(0, 100):
        ind = (i * 5) + 5
        if ind < (l1 - length):
            (s11, s22, p, overlap) = trim(s1, s2, l1, l2, ind, length)
            if p == 1:
                seq = (
                    s1[0 : s1.index(overlap)]
                    + s2[(s2.index(overlap)) : l2].lower()
                )
                if len(seq) > 120:
                    failed = 0
                    break
        # else:break
    return (seq, failed)


def trim(
    s1: str, s2: str, l1: int, l2: int, indent: float, length: int
) -> Tuple[str, str, int, str]:
    """Trim sequences.

    Parameters
    ----------
    s1 : str
        Description
    s2 : str
        Description
    l1 : int
        Description
    l2 : int
        Description
    indent : float
        Description
    length : int
        Description

    Returns
    -------
    Tuple[str, str, int, str]
        trimmed sequences.
    """
    (i, p) = (indent, 0)
    sample1 = s1[i : i + length]
    # print sample1, "\t",s1, s2,"\n"
    index = s2.find(sample1)
    s1a, s2a = s1, s2
    if index != -1:
        if index > i:
            s2a = s2a[index - i : l2]
        else:
            s1a = s1a[i - index : l1]
        min_len = min([len(s1), len(s2)])
        s1a = s1a[0:min_len]
        s2a = s2a[0:min_len]
        p = 1
    return (s1a, s2a, p, sample1)


def check_fasta_not_empty(file: Path) -> int:
    """Checks if fasta is empty.

    Parameters
    ----------
    file : Path
        path to fasta file

    Returns
    -------
    int
        whether file is empty. 1 is not empty, 0 is empty.
    """
    fh = open(file, "r")
    pfh = 0
    for l in fh:
        print(l)
        if len(l) != 0:
            pfh = 1
        break
    fh.close()
    return pfh


def get_freq(header: str) -> int:
    """Get frequncy of isotype assignments.

    Parameters
    ----------
    header : str
        header of sequence.

    Returns
    -------
    int
        frequency
    """
    header = header.split(READ_NUMBER_DIVISION)
    return int(header[1])


def get_match(primer: str, seq: str) -> List:
    """Find sequence in primer.

    Parameters
    ----------
    primer : str
        primer sequence
    seq : str
        sequence

    Returns
    -------
    List
        List of matches.
    """
    loc = []
    if seq.count(primer) != 0:
        for m in re.finditer(primer, seq):
            loc = loc + [m.start()]
            # loc = seq.index(primer)
    return loc


def cluster_i(input_file: Path, clustered_file: Path, identity: float):
    """Run CD-HIT clustering

    Parameters
    ----------
    input_file : Path
        input file to run CD-HIT clustering.
    clustered_file : Path
        output file of CD-HIT.
    identity : float
        sequence identity threshold.
    """
    command = [
        "cd-hit",
        "-i",
        str(input_file),
        "-o",
        str(clustered_file),
        "-c",
        str(identity),
        "-g",
        "1",
        "-d",
        "180",
        "-T",
        "10",
        "-M",
        "0",
        "-AL",
        "40",
        "-bak",
        "1",
        "-p",
        "1",
    ]
    subprocess.run(command)
