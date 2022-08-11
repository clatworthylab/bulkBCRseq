#!/usr/bin/env python
from operator import itemgetter
from pathlib import Path
from typing import Dict, Tuple

from isotyper.utilities import fasta_iterator, cluster_i, write_out, Tree
from isotyper.utilities import READ_NUMBER_DIVISION, EDGE_LENGTHS
from isotyper.qualitycontrol import FILTERED_OUT_NT


def generate_networks(
    sequence_file: Path,
    clustered_file: Path,
    identity: float,
    tmp_file: Path,
    file_out: Path,
):
    """Summary

    Parameters
    ----------
    sequence_file : Path
        input file to run CD-HIT clustering and to retrieve clusters.
    clustered_file : Path
        output file of CD-HIT.
    identity : float
        sequence identity threshold.
    tmp_file : Path
        Description
    file_out : Path
        Description
    """
    cluster_i(
        input_file=sequence_file,
        clustered_file=clustered_file,
        identity=identity,
    )
    s_sizes, cluster = get_cluster_sizes_single(file=clustered_file)
    seqs = get_seqs_single(file=sequence_file)
    get_similar_clusters(s_sizes, cluster, seqs, tmp_file + "_coclustered")
    inv, coclust = get_coclustered(tmp_file + "_coclustered")
    get_cluster_similarities_single(seqs, coclust, cluster, file_out, inv)


def get_clusters(file: Path) -> "Tree":
    """Retrieve clusters from clustered file.

    Parameters
    ----------
    file : Path
        path to clustered file after CD-HIT.

    Returns
    -------
    Tree
        Clustering information in a tree data structure (dictionary).
    """
    fh = open(file, "r")
    cluster = Tree()
    for l in fh:
        l = l.strip()
        l = l.split()
        id = l[2].replace("...", "")
        id = id.replace(">", "")
        cluster[l[0]][id].value = 1
    fh.close()
    return cluster


def get_cluster_sizes_single(file: Path) -> Tuple:
    """Obtain cluster sizes.

    Parameters
    ----------
    file : Path
        path to clustered file

    Returns
    -------
    Tuple
        tuple of size and cluster name
    """
    (cluster) = get_clusters(file=file.with_suffix(file.suffix + ".bak.clstr"))
    print(len(cluster))
    sizes = []
    for c in cluster:
        tot = len(cluster[c])
        t = (c, tot)
        if tot > 50:
            sizes.append(t)
    s = sorted(sizes, key=itemgetter(1), reverse=True)
    return (s, cluster)


def get_seqs_single(file: Path) -> Dict:
    """return dictionary of sequences from fasta file.

    Parameters
    ----------
    file : Path
        fasta file.

    Returns
    -------
    Dict
        sequences as a dictionary.
    """
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split(READ_NUMBER_DIVISION)[0]] = sequence
    fh.close()
    return seqs


def count_diffs(s1: str, s2: str, mis: int) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : str
        Description
    s2 : str
        Description
    mis : int
        Description

    Returns
    -------
    Tuple
        Description
    """
    i1 = 0
    i2 = 0
    mm = 0
    p = 1
    for i in range(0, len(s2) - 1):
        if s1[i1] == s2[i2]:
            i1 = i1 + 1
            i2 = i2 + 1
        else:
            if s1[i1 + 1] == s2[i2]:
                i1 = i1 + 1
                mm = mm + 1
            else:
                if s1[i1] == s2[i2 + 1]:
                    i2 = i2 + 1
                    mm = mm + 1
                else:
                    mm = mm + 1
        if mm > mis:
            p = 0
            break
    return (mm, p)


def get_vaguely_similar_seqs(s1: str, s2: str, mis: int) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : str
        Description
    s2 : str
        Description
    mis : int
        Description

    Returns
    -------
    Tuple
        Description
    """
    l1 = len(s1)
    l2 = len(s2)
    trim = 15
    seg_length = (l1 - (2 * trim)) / mis
    p = 0
    for i in range(0, mis):
        pos = trim + (i * seg_length)
        seg = s1[pos : pos + seg_length]
        index = s2.find(seg)
        if index != -1:
            if index > pos:
                s2 = s2[index - pos : l2]
            else:
                s1 = s1[pos - index : l1]
            min_len = min([len(s1), len(s2)])
            s1 = s1[0:min_len]
            s2 = s2[0:min_len]
            p = 1
            break
    return (s1, s2, p)


def get_similar_clusters(s_sizes, cluster, seqs, tmp_file: Path):
    """Summary

    Parameters
    ----------
    s_sizes : TYPE
        Description
    cluster : TYPE
        Description
    seqs : TYPE
        Description
    tmp_file : TYPE
        Description
    """
    coclust = Tree()
    inv = {}
    mis = 5
    out = ""
    indw = 0
    fh1 = open(tmp_file, "w")
    fh1.close()
    comp = 5
    for i in range(0, len(s_sizes)):
        print(len(s_sizes), i)
        clust = s_sizes[i]
        c1 = clust[0]
        if c1 not in inv:
            s1 = ""
            coclust[c1][c1].value = 1
            inv[c1] = c1
            seqs1 = []
            ind = 0
            for seq_id in cluster[c1]:
                if seq_id in seqs:
                    s1 = seqs[seq_id]
                    seqs1.append(s1)
                    ind = ind + 1
                    if ind > comp:
                        break
            for i2 in range(i + 1, len(s_sizes)):
                clust2 = s_sizes[i2]
                c2 = clust2[0]
                ind = 0
                found = 0
                for seq_id in cluster[c2]:
                    if seq_id in seqs:
                        # Get difference with c1 sequences
                        # If sequence is different > X% then move onto next cluster
                        # If sequence is less different than Y% then co-cluster the clusters
                        p = 0
                        ind = ind + 1
                        for s1 in seqs1:
                            (s1, s2, p) = get_vaguely_similar_seqs(
                                s1, seqs[seq_id], mis
                            )
                            if p == 1:
                                (mm, q) = count_diffs(s1, s2, mis)
                                if q == 1:
                                    out = out + c1 + "\t" + c2 + "\n"
                                    indw = indw + 1
                                    if indw > 100:
                                        write_out(out, tmp_file)
                                        out = ""
                                        indw = 0
                                    found = 1
                                    break
                        if ind > comp or found == 1:
                            break
    write_out(out, tmp_file)
    print(out)


def main():
    """main function in step 3."""
    generate_networks(
        sequence_file=FILTERED_OUT_NT,
        clustered_file=CLUSTERED_TMP_FILE,
        identity=EDGE_LENGTHS,
        tmp_file=tmp_file,
        file_out=att_file,
    )


if __name__ == "__main__":
    main()
