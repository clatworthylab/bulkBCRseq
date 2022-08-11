#!/usr/bin/env python
import networkx as nx

from operator import itemgetter
from pathlib import Path
from typing import Dict, Tuple, List

from isotyper.utilities import (
    create_file,
    cluster_i,
    fasta_iterator,
    Tree,
    write_out,
    EDGE_LENGTHS,
    READ_NUMBER_DIVISION,
)
from isotyper.qualitycontrol import FILTERED_OUT_NT
from isotyper.network._settings import (
    ATT_FILE,
    CHKEDGES_FILE,
    CLUSTERED_TMP_FILE,
    CLUSTER_FILE,
    COCLUSTERED,
    EDGES_FILE,
    PLOT_IDS_FILE,
    REDUCED_FILE,
    SEQS_FILE,
    VERTEX_REL_FILE,
)
from isotyper.network._decon import deconvolute_edges


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
            i1 += 1
            i2 += 1
        else:
            if s1[i1 + 1] == s2[i2]:
                i1 += 1
                mm += 1
            else:
                if s1[i1] == s2[i2 + 1]:
                    i2 += 1
                    mm += 1
                else:
                    mm += 1
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
        sequence 1
    s2 : str
        sequence 2
    mis : int
        mismatch

    Returns
    -------
    Tuple[str, str, int]
        vaguely similar sequences.
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
                    ind += 1
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
                        ind += 1
                        for s1 in seqs1:
                            (s1, s2, p) = get_vaguely_similar_seqs(
                                s1=s1, s2=seqs[seq_id], mis=mis
                            )
                            if p == 1:
                                (mm, q) = count_diffs(s1=s1, s2=s2, mis=mis)
                                if q == 1:
                                    out = out + c1 + "\t" + c2 + "\n"
                                    indw += 1
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
            i1 += 1
            i2 += 1
        else:
            mm += 1
            i2 += 1
            i1 += 1
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
                            ind1 += 1
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
                                ind1 += 1
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
    create_file(file_out)
    ind = 0
    total = len(cluster)
    t = 0
    for c in cluster:
        if c not in inv:
            ind += 1
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
                t += 1
            for c1 in coclust[c]:
                for id in cluster[c1]:
                    clust_seqs.append(
                        (
                            id,
                            seqs[id.split(READ_NUMBER_DIVISION)[0]],
                            len(seqs[id.split(READ_NUMBER_DIVISION)[0]]),
                        )
                    )
                    t += 1
            clust_seqs = sorted(clust_seqs, key=itemgetter(2), reverse=True)
            if len(clust_seqs) > 1:
                if len(clust_seqs) > 500:
                    print(len(clust_seqs), total, ind)
                get_similarity_single(
                    clust_seqs=clust_seqs, file_out=file_out, c=c
                )


def read_graphical_inputs(
    att_file: Path, file_edges: Path
) -> Tuple[nx.Graph, int]:
    """Load network into networkx.

    Parameters
    ----------
    att_file : Path
        path to Att.txt (final node table)
    file_edges : Path
        path to Edges.txt (final edge table)

    Returns
    -------
    Tuple[nx.Graph, int]
        networkx graph holding the BCR network.
    """
    fh = open(att_file, "r")
    freq = {}
    size = []
    ind = 0
    G = nx.Graph()
    G.rtt = {}
    for l in fh:
        l = l.strip()
        l = l.split()
        f = int(l[1])
        freq[l[0]] = f
        size.append(f)
        G.add_node(l[0])
        G.rtt[l[0]] = int(f)
        ind += 1
    scale1 = max(size)
    fh.close()
    fh = open(file_edges, "r")
    for l in fh:
        l = l.strip()
        l = l.split()
        if l[0] in freq and l[1] in freq:
            G.add_edge(l[0], l[1])
    fh.close()
    return (G, scale1)


def output_cluster_file(graph: nx.Graph, cluster_file: Path):
    """Write out cluster file.

    Parameters
    ----------
    graph : nx.Graph
        networkx graph.
    cluster_file : Path
        path to Cluster_sample_identities.txt
    """
    con = nx.connected_components(graph)
    ind = 0
    ind1 = 0
    ind2 = 0
    create_file(cluster_file)
    out = "# Connected_components\n"
    max_f, t, nvertmax = 0, 0, 0
    for i in con:
        ind += 1
        tc = 0
        nvert = 0
        for j in i:
            ind1 += 1
            ind2 += 1
            out = (
                out
                + str(ind1)
                + "\t"
                + str(ind)
                + "\t"
                + j
                + "\t"
                + str(graph.rtt[j])
                + "\n"
            )
            tc, t = tc + graph.rtt[j], t + graph.rtt[j]
            nvert += 1
            if ind2 > 100:
                write_out(out, cluster_file)
                out = ""
                ind2 = 0
        if tc > max_f:
            max_f, nvertmax = tc, nvert
    write_out(out, cluster_file)
    print(
        file_vertex,
        "Maximum cluster:",
        max_f * 100.0 / t,
        "%",
        nvertmax,
        "vertices",
    )


def get_network_clusters(att_file: Path, file_edges: Path, cluster_file: Path):
    """Get network clusters.

    Parameters
    ----------
    att_file : Path
        path to Att.txt (final node table)
    file_edges : Path
        path to Edges.txt (final edge table)
    cluster_file : Path
        path to Cluster_sample_identities.txt
    """
    (G, scale) = read_graphical_inputs(att_file=att_file, file_edges=file_edges)
    output_cluster_file(graph=G, cluster_file=cluster_file)


def reduce_identical_sequences(reduced_file: Path, att_file: Path):
    """Summary

    Parameters
    ----------
    reduced_file : Path
        path to Fully reduced fasta file.
    att_file : Path
        path to Att.txt (final node table)
    """
    create_file(reduced_file)
    fh = open(att_file, "r")
    (ind, out) = (0, "")
    for l in fh:
        l = l.strip().split()
        if l[0].count("|") >= 1:
            out = out + ">" + l[0] + "\n" + l[2] + "\n"
        else:
            out = (
                out
                + ">"
                + l[0].split(READ_NUMBER_DIVISION)[0]
                + READ_NUMBER_DIVISION
                + l[1]
                + "\n"
                + l[2]
                + "\n"
            )
        ind += 1
        if ind > 500:
            write_out(out, reduced_file)
            (ind, out) = (0, "")
    fh.close()
    write_out(out, reduced_file)


def get_network_input(
    file_cluster: Path, outfile: Path, edge_file: Path, checked_edges: Path
):
    """Summary

    Parameters
    ----------
    file_cluster : Path
        Description
    outfile : Path
        path to Plot_ids.txt
    edge_file : Path
        path to Edges.txt (final edge table)
    checked_edges : Path
        path to Checked_edges.txt
    """
    create_file(outfile)
    fh = open(file_cluster, "r")
    cluster = Tree()
    ids = {}
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            cluster[l[1]][l[1] + "\t" + l[2] + "\t" + l[3]].value = 1
            ids[l[2]] = 1
    fh.close()
    (out, ind) = ("", 0)
    for c in cluster:
        if len(cluster[c]) > 1:
            for l in cluster[c]:
                out = out + l + "\n"
                ind += 1
                if int(l.split("\t")[2]) > 1000:
                    print(l)
        else:
            for l in cluster[c]:
                if int(l.split("\t")[2]) > 1:
                    out = out + l + "\n"
                    ind += 1
        if ind > 300:
            write_out(out, outfile)
            (out, ind) = ("", 0)
    write_out(out, outfile)
    create_file(checked_edges)
    fh = open(edge_file, "r")
    (out, ind) = ("", 0)
    for l in fh:
        l = l.strip().split()
        if l[0] in ids and l[1] in ids:
            out = out + l[0] + "\t" + l[1] + "\t" + l[2] + "\n"
            ind += 1
            if ind > 300:
                write_out(out, checked_edges)
                (out, ind) = ("", 0)
    write_out(out, checked_edges)
    fh.close()


def main():
    """main function in step 3."""
    generate_networks(
        sequence_file=FILTERED_OUT_NT,
        clustered_file=CLUSTERED_TMP_FILE,
        identity=EDGE_LENGTHS,
        coclustered_file=COCLUSTERED,
        file_out=VERTEX_REL_FILE,
    )
    deconvolute_edges(
        seq_file=FILTERED_OUT_NT,
        file_vertex=VERTEX_REL_FILE,
        att_file=ATT_FILE,
        file_seqs=SEQS_FILE,
        file_edges=EDGES_FILE,
    )
    get_network_clusters(
        att_file=ATT_FILE, file_edge=EDGES_FILE, cluster_file=CLUSTER_FILE
    )
    reduce_identical_sequences(reduced_file=REDUCED_FILE, att_file=ATT_FILE)
    get_network_input(
        file_cluster=CLUSTER_FILE,
        outfile=PLOT_IDS_FILE,
        edge_file=EDGES_FILE,
        checked_edges=CHKEDGES_FILE,
    )


if __name__ == "__main__":
    main()
