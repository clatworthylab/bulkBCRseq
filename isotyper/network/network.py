#!/usr/bin/env python
from operator import itemgetter
from pathlib import Path
from typing import Dict, Tuple, List

from isotyper.utilities import fasta_iterator, cluster_i, write_out, Tree
from isotyper.utilities import READ_NUMBER_DIVISION, EDGE_LENGTHS
from isotyper.qualitycontrol import FILTERED_OUT_NT
from isotyper.network._settings import COCLUSTERED, CLUSTERED_TMP_FILE, ATT_FILE


def generate_networks(
    sequence_file: Path,
    clustered_file: Path,
    identity: float,
    coclustered_file: Path,
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
    coclustered_file : Path
        path to coclustered file.
    file_out : Path
        out path for vertices table.
    """
    cluster_i(
        input_file=sequence_file,
        clustered_file=clustered_file,
        identity=identity,
    )
    s_sizes, cluster = get_cluster_sizes_single(file=clustered_file)
    seqs = get_seqs_single(file=sequence_file)
    get_similar_clusters(
        s_sizes=s_sizes,
        cluster=cluster,
        seqs=seqs,
        coclustered_file=coclustered_file,
    )
    inv, coclust = get_coclustered(file=coclustered_file)
    get_cluster_similarities_single(
        seqs=seqs, coclust=coclust, cluster=cluster, file_out=file_out, inv=inv
    )


def get_clusters(file: Path) -> Tree:
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


def get_cluster_sizes_single(file: Path) -> Tuple[List, Tree]:
    """Obtain cluster sizes.

    Parameters
    ----------
    file : Path
        path to clustered file

    Returns
    -------
    Tuple[List, Tree]
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


def count_diffs(s1: str, s2: str, mis: int) -> Tuple[int, int]:
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
    Tuple[int, int]
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


def get_vaguely_similar_seqs(
    s1: str, s2: str, mis: int
) -> Tuple[str, str, int]:
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
    Tuple[str, str, int]
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


def get_similar_clusters(
    s_sizes: List, cluster: Tree, seqs: Dict, coclustered_file: Path
):
    """Get similar clusters.

    Parameters
    ----------
    s_sizes : List
        size of clusters.
    cluster : Tree
        clusters of sequences.
    seqs : Dict
        dictionary of sequences.
    coclustered_file : Path
        output coclustered file.
    """
    coclust = Tree()
    inv = {}
    mis = 5
    out = ""
    indw = 0
    fh1 = open(coclustered_file, "w")
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
                                        write_out(out, coclustered_file)
                                        out = ""
                                        indw = 0
                                    found = 1
                                    break
                        if ind > comp or found == 1:
                            break
    write_out(out, coclustered_file)
    print(out)


def get_coclustered(file: Path) -> Tuple[Dict, Tree]:
    """retrieve coclustered sequences and clusters.

    Parameters
    ----------
    file : Path
        path to coclustered file.

    Returns
    -------
    Tuple[Dict, Tree]
        dictionary of clustered sequences
    """
    inv = {}
    coclust = Tree()
    fh = open(file, "r")
    for l in fh:
        l = l.strip()
        l = l.split()
        inv[l[1]] = l[0]
        coclust[l[0]][l[1]].value = 1
    fh.close()
    return (inv, coclust)


def trim_sequences(s1: str, s2: str, l1: int, l2: int) -> Tuple[str, str, int]:
    """trim sequences

    Parameters
    ----------
    s1 : str
        sequence 1
    s2 : str
        sequence 2
    l1 : int
        length 1
    l2 : int
        length 2

    Returns
    -------
    Tuple[str, str, int]
        trimmed sequences.
    """
    p = 0
    i = 15
    sample1 = s1[i : i + 25]
    index = s2.find(sample1)
    if index != -1:
        if index > i:
            if index - i <= 20:
                s2 = s2[index - i : l2]
        else:
            if i - index <= 20:
                s1 = s1[i - index : l1]
        min_len = min([len(s1), len(s2)])
        if (max([len(s1), len(s2)]) - min_len) < 25:
            s1 = s1[0:min_len]
            s2 = s2[0:min_len]
            p = 1
    else:
        i = l1 - 50
        sample1 = s1[i : i + 25]
        index = s2.find(sample1)
        if index != -1:
            if index > i:
                if index - i <= 20:
                    s2 = s2[index - i : l2]
            else:
                if i - index <= 20:
                    s1 = s1[i - index : l1]
            min_len = min([len(s1), len(s2)])
            if (max([len(s1), len(s2)]) - min_len) < 25:
                s1 = s1[0:min_len]
                s2 = s2[0:min_len]
                p = 1
            else:
                p = 0
        else:
            p = 0
    return (s1, s2, p)


def do_counting(s1: str, s2: str, mismatch: int) -> int:
    """Summary

    Parameters
    ----------
    s1 : str
        sequence 1
    s2 : str
        sequence 2
    mismatch : int
        mismatch threshold

    Returns
    -------
    int
        mismatch size.
    """
    i1 = 0
    i2 = 0
    mm = 0
    l = min([len(s1), len(s2)])
    for i in range(0, l - 2):
        if s1[i1] == s2[i2]:
            i1 = i1 + 1
            i2 = i2 + 1
        else:
            mm = mm + 1
            i2 = i2 + 1
            i1 = i1 + 1
        if mm > mismatch + 1:
            break
    return mm


def get_diff(s1: str, s2: str, mismatch: int) -> Tuple[int, int]:
    """Summary

    Parameters
    ----------
    s1 : str
        sequence 1
    s2 : str
        sequence 2
    mismatch : int
        mismatch size

    Returns
    -------
    Tuple[int, int]
        boolean, mismatch
    """
    p = 1
    mm = do_counting(s1, s2, mismatch)
    if mm > mismatch:
        p = 0
    return (p, mm)


def get_similarity_single(clust_seqs: List, file_out: Path, c: str):
    """get similarity in clusters.

    Parameters
    ----------
    clust_seqs : List
        clusters
    file_out : Path
        path to output file.
    c : str
        cluster name
    """
    done = Tree()
    total = len(clust_seqs)
    out = ""
    ind1 = 0
    mismatch = 1
    for i1 in range(0, total - 1):
        read1 = clust_seqs[i1]
        id1 = read1[0]
        seq1 = read1[1]
        l1 = read1[2]
        done[id1] = 1
        for i2 in range(i1 + 1, total):
            if i1 < i2:
                read2 = clust_seqs[i2]
                id2 = read2[0]
                seq2 = read2[1]
                l2 = read2[2]
                if id2 not in done:
                    (s1, s2, p1) = trim_sequences(
                        s1=seq1, s2=seq2, l1=l1, l2=l2
                    )
                    if p1 == 1:
                        if s1 == s2:
                            out = (
                                out
                                + "0\t"
                                + id1
                                + "\t"
                                + id2
                                + "\t"
                                + c
                                + "\t"
                                + str(l1)
                                + "\t"
                                + str(l2)
                                + "\n"
                            )
                            ind1 = ind1 + 1
                        else:
                            (p, mm) = get_diff(s1=s1, s2=s2, mismatch=mismatch)
                            if p == 1 and mm <= mismatch:
                                out = (
                                    out
                                    + str(mm)
                                    + "\t"
                                    + id1
                                    + "\t"
                                    + id2
                                    + "\t"
                                    + c
                                    + "\t"
                                    + str(l1)
                                    + "\t"
                                    + str(l2)
                                    + "\n"
                                )
                                ind1 = ind1 + 1
                if ind1 >= 100:
                    write_out(out, file_out)
                    ind1 = 0
                    out = ""
    write_out(out, file_out)


def get_cluster_similarities_single(
    seqs: Dict, coclust: Tree, cluster: Tree, file_out: Path, inv: Dict
):
    """Summary

    Parameters
    ----------
    seqs : Dict
        sequences
    coclust : Tree
        coclustered sequences
    cluster : Tree
        clusters
    file_out : Path
        out path
    inv : Dict
        dictionary
    """
    fh = open(file_out, "w")
    fh.close()
    ind = 0
    total = len(cluster)
    t = 0
    for c in cluster:
        if c not in inv:
            ind = ind + 1
            clust_seqs = []
            t = 0
            for id in cluster[c]:
                clust_seqs.append(
                    (
                        id,
                        seqs[id.split(READ_NUMBER_DIVISION)[0]],
                        len(seqs[id.split(READ_NUMBER_DIVISION)[0]]),
                    )
                )
                t = t + 1
            for c1 in coclust[c]:
                for id in cluster[c1]:
                    clust_seqs.append(
                        (
                            id,
                            seqs[id.split(READ_NUMBER_DIVISION)[0]],
                            len(seqs[id.split(READ_NUMBER_DIVISION)[0]]),
                        )
                    )
                    t = t + 1
            clust_seqs = sorted(clust_seqs, key=itemgetter(2), reverse=True)
            if len(clust_seqs) > 1:
                if len(clust_seqs) > 500:
                    print(len(clust_seqs), total, ind)
                get_similarity_single(
                    clust_seqs=clust_seqs, file_out=file_out, c=c
                )


def main():
    """main function in step 3."""
    generate_networks(
        sequence_file=FILTERED_OUT_NT,
        clustered_file=CLUSTERED_TMP_FILE,
        identity=EDGE_LENGTHS,
        coclustered_file=COCLUSTERED,
        file_out=ATT_FILE,
    )
    deconvolute_edges(
        seq_file=FILTERED_OUT_NT,
        att_file=ATT_FILE,
        file_vertex=file_vertex,
        file_seqs=file_seqs,
        tmp_file0=tmp_file0,
        file_edges=file_edges,
        read_number_division=READ_NUMBER_DIVISION,
    )


if __name__ == "__main__":
    main()
