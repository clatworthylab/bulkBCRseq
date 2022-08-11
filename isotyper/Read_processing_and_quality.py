#!/usr/bin/env python
import os
import re
import shutil
import subprocess
import sys

import networkx as nx

from glob import glob
from operator import itemgetter, add
from pathlib import Path
from typing import Union, Tuple, Dict

from isotyper.utilities import *
from isotyper.qualitycontrol import *
from isotyper.network import *


def get_similarity_single(clust_seqs, file_out, c):
    """Summary

    Parameters
    ----------
    clust_seqs : TYPE
        Description
    file_out : TYPE
        Description
    c : TYPE
        Description
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
                    (s1, s2, p1) = trim_sequences(seq1, seq2, l1, l2)
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
                            (p, mm) = get_diff(s1, s2, mismatch)
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


def get_diff(s1, s2, mismatch) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    mismatch : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    p = 1
    (mm) = do_counting(s1, s2, mismatch)
    if mm > mismatch:
        p = 0
    return (p, mm)


def do_counting(s1, s2, mismatch):
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    mismatch : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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


def trim_sequences(s1, s2, l1, l2) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    l1 : TYPE
        Description
    l2 : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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


def get_cluster_similarities_single(seqs, coclust, cluster, file_out, inv):
    """Summary

    Parameters
    ----------
    seqs : TYPE
        Description
    coclust : TYPE
        Description
    cluster : TYPE
        Description
    file_out : TYPE
        Description
    inv : TYPE
        Description
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
                get_similarity_single(clust_seqs, file_out, c)


def get_clusters(file: Path) -> "Tree":
    """Retrieve clusters from file.

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
    """Summary

    Parameters
    ----------
    file : Path
        Description

    Returns
    -------
    Tuple
        Description
    """
    (cluster) = get_clusters(file + ".bak.clstr")
    print(len(cluster))
    sizes = []
    for c in cluster:
        tot = len(cluster[c])
        t = (c, tot)
        if tot > 50:
            sizes.append(t)
    s = sorted(sizes, key=itemgetter(1), reverse=True)
    return (s, cluster)


def get_vaguely_similar_seqs(s1, s2, mis) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    mis : TYPE
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


def count_diffs(s1, s2, mis) -> Tuple:
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    mis : TYPE
        Description

    Returns
    -------
    TYPE
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


def get_similar_clusters(s_sizes, cluster, seqs, tmp_file):
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
            for id in cluster[c1]:
                if id in seqs:
                    s1 = seqs[id]
                    seqs1.append(s1)
                    ind = ind + 1
                    if ind > comp:
                        break
            for i2 in range(i + 1, len(s_sizes)):
                clust2 = s_sizes[i2]
                c2 = clust2[0]
                ind = 0
                found = 0
                for id in cluster[c2]:
                    if id in seqs:
                        # Get difference with c1 sequences
                        # If sequence is different > X% then move onto next cluster
                        # If sequence is less different than Y% then co-cluster the clusters
                        p = 0
                        ind = ind + 1
                        for s1 in seqs1:
                            (s1, s2, p) = get_vaguely_similar_seqs(
                                s1, seqs[id], mis
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


def get_coclustered(file) -> Tuple:
    """Summary

    Parameters
    ----------
    file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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


def decon_edges(att_file, file_seqs, file_edges, file_vertex, seq_file):
    """Summary

    Parameters
    ----------
    att_file : TYPE
        Description
    file_seqs : TYPE
        Description
    file_edges : TYPE
        Description
    file_vertex : TYPE
        Description
    seq_file : TYPE
        Description
    """
    fh = open(seq_file, "r")
    freq_id = {}
    for header, sequence in fasta_iterator(fh):
        freq_id[header.split(READ_NUMBER_DIVISION)[0]] = header
    fh.close()
    inverse, raw = get_inverse_ids(file_seqs, file_vertex, freq_id)
    edges, edges23 = Tree(), Tree()
    fh1 = open(att_file, "r")
    for l in fh1:
        l = l.strip()
        l1 = l.split()
        if int(l1[0]) == 1 or int(l1[0]) == 2:
            id1, id2 = (
                l1[1].split(READ_NUMBER_DIVISION)[0],
                l1[2].split(READ_NUMBER_DIVISION)[0],
            )
            id1, id2 = freq_id[id1], freq_id[id2]
            edges[id2][id1].value = 1
    fh1.close()
    print_single_edges(file_edges, inverse, edges, tmp_file_1, raw)


def get_inverse_ids(file_seqs, file_vertex, freq_id):
    """Summary

    Parameters
    ----------
    file_seqs : TYPE
        Description
    file_vertex : TYPE
        Description
    freq_id : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(file_seqs, "r")
    inverse, raw = {}, {}
    ind = 0
    for l in fh:
        if l[0] == ">":
            l = l.strip()
            l = l.replace(">", "")
            l = l.split("||")
            inverse[freq_id[l[0].split(READ_NUMBER_DIVISION)[0]]] = freq_id[
                l[1].split(READ_NUMBER_DIVISION)[0]
            ]
            raw[l[1].split(READ_NUMBER_DIVISION)[0]] = 1
            ind = ind + 1
    fh.close()
    fh = open(file_vertex, "r")
    for l in fh:
        l = l.strip().split()
        if l[0].split(READ_NUMBER_DIVISION)[0] in raw:
            raw[l[0].split(READ_NUMBER_DIVISION)[0]] = l[0]
    fh.close()
    return (inverse, raw)


def decon_identical(
    seq_file, att_file, file_vertex, file_seqs, tmp_file, read_number_division
):
    """Summary

    Parameters
    ----------
    seq_file : TYPE
        Description
    att_file : TYPE
        Description
    file_vertex : TYPE
        Description
    file_seqs : TYPE
        Description
    tmp_file : TYPE
        Description
    read_number_division : TYPE
        Description
    """
    fh = open(seq_file, "r")
    all, seqs, freq_id = {}, {}, {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split(READ_NUMBER_DIVISION)[0]] = sequence
        all[header.split(READ_NUMBER_DIVISION)[0]] = sequence
        freq_id[header.split(READ_NUMBER_DIVISION)[0]] = header
    fh.close()
    fh1 = open(att_file, "r")
    same1 = Tree()
    for l in fh1:
        l = l.strip()
        l1 = l.split()
        cluster = l1[3]
        if l1[0] == "0":
            same1[cluster][l1[2]][l1[1]].value = 1
            same1[cluster][l1[1]][l1[2]].value = 1
    fh1.close()
    same, inverse, out, ind, length = {}, {}, "", 0, {}
    fh = open(file_seqs, "w")
    fh1.close()
    j = header
    for c in same1:
        sub_same = Tree()
        for id1 in same1[c]:
            for id2 in same1[c][id1]:
                sub_same[id1.split(READ_NUMBER_DIVISION)[0]][
                    id2.split(READ_NUMBER_DIVISION)[0]
                ].value = 1
                sub_same[id2.split(READ_NUMBER_DIVISION)[0]][
                    id1.split(READ_NUMBER_DIVISION)[0]
                ].value = 1
        (sub_same, sub_inv) = deconvolute_same_array(sub_same)
        for i in sub_same:
            s = seqs[i.split(READ_NUMBER_DIVISION)[0]]
            total = 0
            mins = s
            for j in sub_same[i]:
                j1 = freq_id[j.split(READ_NUMBER_DIVISION)[0]]
                freq = list(
                    map(
                        int,
                        j1.split(read_number_division)[1]
                        .split("|")[0]
                        .split("_"),
                    )
                )
                if total == 0:
                    total = freq
                else:
                    total = list(map(add, freq, total))
                inverse[j.split(READ_NUMBER_DIVISION)[0]] = i
                out = (
                    out
                    + ">"
                    + j
                    + "||"
                    + i
                    + "\n"
                    + seqs[j.split(READ_NUMBER_DIVISION)[0]]
                    + "\n"
                )
                s = seqs[j.split(READ_NUMBER_DIVISION)[0]]
                if len(s) < len(mins):
                    mins = s
                ind = ind + 1
                if ind > 100:
                    write_out(out, file_seqs)
                    out = ""
                    ind = 0
            same[i] = total
            length[i] = mins
    info = ""
    if (
        len(
            freq_id[j.split(READ_NUMBER_DIVISION)[0]]
            .split(read_number_division)[1]
            .split("|")
        )
        >= 2
    ):
        info = (
            "|"
            + freq_id[j.split(READ_NUMBER_DIVISION)[0]]
            .split(read_number_division)[1]
            .split("|")[1]
        )
    write_out(out, file_seqs)
    print_vertices(
        all,
        inverse,
        same,
        file_vertex,
        length,
        read_number_division,
        info,
        freq_id,
    )


def print_single_edges(file_edges, inverse, edges, tmp_file, raw):
    """Summary

    Parameters
    ----------
    file_edges : TYPE
        Description
    inverse : TYPE
        Description
    edges : TYPE
        Description
    tmp_file : TYPE
        Description
    raw : TYPE
        Description
    """
    fh = open(tmp_file, "w")
    fh.close()
    edge, ind = "", 0
    for id1 in edges:
        ida = id1
        if id1 in inverse:
            ida = inverse[id1]
            if ida.split(READ_NUMBER_DIVISION)[0] in raw:
                ida = raw[ida.split(READ_NUMBER_DIVISION)[0]]
        for id2 in edges[id1]:
            idb = id2
            if id2 in inverse:
                idb = inverse[id2]
                if idb.split(READ_NUMBER_DIVISION)[0] in raw:
                    idb = raw[idb.split(READ_NUMBER_DIVISION)[0]]
            if ida != idb:
                edge = (
                    edge
                    + ida
                    + "\t"
                    + idb
                    + "\t"
                    + str(1)
                    + "\t"
                    + id1
                    + "\t"
                    + id2
                    + "\n"
                )
                ind = ind + 1
                if ind > 300:
                    write_out(edge, tmp_file)
                    edge = ""
                    ind = 0
    write_out(edge, tmp_file)


def print_vertices(
    all, inverse, same, file_vertex, length, read_number_division, info, freq_id
):
    """Summary

    Parameters
    ----------
    all : TYPE
        Description
    inverse : TYPE
        Description
    same : TYPE
        Description
    file_vertex : TYPE
        Description
    length : TYPE
        Description
    read_number_division : TYPE
        Description
    info : TYPE
        Description
    freq_id : TYPE
        Description
    """
    out = ""
    # total = 0
    fh = open(file_vertex, "w")
    fh.close()
    ind = 0
    for id in all:
        ind = ind + 1
        if id in inverse:
            if id in same:
                id1 = (
                    id.split(read_number_division)[0]
                    + read_number_division
                    + "_".join(map(str, same[id]))
                    + info
                )
                if id not in length:
                    out = (
                        out + id1 + "\t" + str(same[id]) + "\t" + all[id] + "\n"
                    )
                else:
                    out = (
                        out
                        + id1
                        + "\t"
                        + str(sum(same[id]))
                        + "\t"
                        + length[id]
                        + "\n"
                    )
        else:
            if id in freq_id:
                id1 = freq_id[id]
            else:
                print(id)
            freq = sum(
                map(
                    int,
                    id1.split(read_number_division)[1].split("|")[0].split("_"),
                )
            )
            out = out + id1 + "\t" + str(freq) + "\t" + all[id] + "\n"
        if ind > 300:
            write_out(out, file_vertex)
            ind = 0
            out = ""
    write_out(out, file_vertex)


def reduce_edges(file_in, file_out):
    """Summary

    Parameters
    ----------
    file_in : TYPE
        Description
    file_out : TYPE
        Description
    """
    done, out, ind = Tree(), "", 0
    fh = open(file_in, "r")
    fh1 = open(file_out, "w")
    fh1.close()
    for l in fh:
        l = l.strip()
        l = l.split()
        if l[0] not in done[l[1]] and l[1] not in done[l[0]]:
            if l[0] != l[1]:
                out = out + l[0] + "\t" + l[1] + "\t" + l[2] + "\n"
                done[l[0]][l[1]].value = 1
                ind = ind + 1
                if ind > 300:
                    ind = 0
                    write_out(out, file_out)
                    out = ""
    write_out(out, file_out)


def deconvolute_edges(
    seq_file,
    att_file,
    file_vertex,
    file_seqs,
    tmp_file0,
    file_edges,
    read_number_division,
):
    """Summary

    Parameters
    ----------
    seq_file : TYPE
        Description
    att_file : TYPE
        Description
    file_vertex : TYPE
        Description
    file_seqs : TYPE
        Description
    tmp_file0 : TYPE
        Description
    file_edges : TYPE
        Description
    read_number_division : TYPE
        Description
    """
    decon_identical(
        seq_file,
        att_file,
        file_vertex,
        file_seqs,
        tmp_file0,
        read_number_division,
    )
    decon_edges(att_file, file_seqs, file_edges, file_vertex, seq_file)
    reduce_edges(tmp_file_1, file_edges)


def get_network_clusters(file_vertex, file_edges, cluster_file):
    """Summary

    Parameters
    ----------
    file_vertex : TYPE
        Description
    file_edges : TYPE
        Description
    cluster_file : TYPE
        Description
    """
    (G, scale) = read_graphical_inputs(file_vertex, file_edges)
    output_cluster_file(G, cluster_file)


def read_graphical_inputs(file_vertex, file_edges):
    """Summary

    Parameters
    ----------
    file_vertex : TYPE
        Description
    file_edges : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(file_vertex, "r")
    freq = {}
    size = []
    ind = 0
    G = nx.Graph()
    G.rtt = {}
    for l in fh:
        l = l.strip()
        l = l.split()
        f = int(
            l[1]
        )  # sum(map(int, l[1].split(read_number_division)[1].split("|")[0].split("_")))
        freq[l[0]] = f
        size.append(f)
        G.add_node(l[0])
        G.rtt[l[0]] = int(f)
        ind = ind + 1
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


def output_cluster_file(G, cluster_file):
    """Summary

    Parameters
    ----------
    G : TYPE
        Description
    cluster_file : TYPE
        Description
    """
    con = nx.connected_components(G)
    ind = 0
    ind1 = 0
    ind2 = 0
    fh = open(cluster_file, "w")
    fh.close()
    out = "# Connected_components\n"
    max_f, t, nvertmax = 0, 0, 0
    for i in con:
        ind = ind + 1
        tc = 0
        nvert = 0
        for j in i:
            ind1 = ind1 + 1
            ind2 = ind2 + 1
            out = (
                out
                + str(ind1)
                + "\t"
                + str(ind)
                + "\t"
                + j
                + "\t"
                + str(G.rtt[j])
                + "\n"
            )
            tc, t = tc + G.rtt[j], t + G.rtt[j]
            nvert = nvert + 1
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


def get_network_input(file_cluster, outfile, edge_file, checked_edges):
    """Summary

    Parameters
    ----------
    file_cluster : TYPE
        Description
    outfile : TYPE
        Description
    edge_file : TYPE
        Description
    checked_edges : TYPE
        Description
    """
    fh = open(outfile, "w")
    fh.close()
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
                ind = ind + 1
                if int(l.split("\t")[2]) > 1000:
                    print(l)
        else:
            for l in cluster[c]:
                if int(l.split("\t")[2]) > 1:
                    out = out + l + "\n"
                    ind = ind + 1
        if ind > 300:
            write_out(out, outfile)
            (out, ind) = ("", 0)
    write_out(out, outfile)
    fh = open(checked_edges, "w")
    fh.close()
    fh = open(edge_file, "r")
    (out, ind) = ("", 0)
    for l in fh:
        l = l.strip().split()
        if l[0] in ids and l[1] in ids:
            out = out + l[0] + "\t" + l[1] + "\t" + l[2] + "\n"
            ind = ind + 1
            if ind > 300:
                write_out(out, checked_edges)
                (out, ind) = ("", 0)
    write_out(out, checked_edges)
    fh.close()


def reduce_identical_sequences(Reduced_file, file_vertex, read_number_division):
    """Summary

    Parameters
    ----------
    Reduced_file : TYPE
        Description
    file_vertex : TYPE
        Description
    read_number_division : TYPE
        Description
    """
    fh = open(Reduced_file, "w")
    fh.close()
    fh = open(file_vertex, "r")
    (ind, out) = (0, "")
    for l in fh:
        l = l.strip().split()
        if l[0].count("|") >= 1:
            out = out + ">" + l[0] + "\n" + l[2] + "\n"
        else:
            out = (
                out
                + ">"
                + l[0].split(read_number_division)[0]
                + read_number_division
                + l[1]
                + "\n"
                + l[2]
                + "\n"
            )
        ind = ind + 1
        if ind > 500:
            write_out(out, Reduced_file)
            (ind, out) = (0, "")
    fh.close()
    write_out(out, Reduced_file)


def generate_networks(
    Sequence_file, tmp_file_1, edge_lengths, tmp_file, file_out
):
    """Summary

    Parameters
    ----------
    Sequence_file : TYPE
        Description
    tmp_file_1 : TYPE
        Description
    edge_lengths : TYPE
        Description
    tmp_file : TYPE
        Description
    file_out : TYPE
        Description
    """
    cluster_i(Sequence_file, tmp_file_1, edge_lengths)
    (s_sizes, cluster) = get_cluster_sizes_single(tmp_file_1)
    (seqs) = get_seqs_single(Sequence_file)
    get_similar_clusters(s_sizes, cluster, seqs, tmp_file + "_coclustered")
    (inv, coclust) = get_coclustered(tmp_file + "_coclustered")
    get_cluster_similarities_single(seqs, coclust, cluster, file_out, inv)


###########################
from isotyper.utilities._args import (
    ORG,
    OUT_FASTQ,
    OUT_NET,
    OUT_ORTSEQ,
    OUT_ORTSEQ_TMP,
    OUT_PATH,
    REVERSE_PRIMER_GROUP,
    SAMPLE_ID,
    SOURCE,
)

# out_path = Path(sys.argv[1])
# sample_id = sys.argv[2]
# barcode_group = sys.argv[3]
# gene = sys.argv[4]
# paired = "OVERLAPPING"
# species = sys.argv[6]
# source = sys.argv[7]
# length = sys.argv[8]
# primer_file = sys.argv[9]
# method = sys.argv[10]
command_source = sys.argv[11]
command_source = command_source.split(",")
# if len(sys.argv) > 13:
#     reverse_primer_group = sys.argv[13]
# else:
#     reverse_primer_group = "OTHER"
print("Reverse primer group: ", REVERSE_PRIMER_GROUP)

# Files for QC and filtering
# seq_file1 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_1.fasta"
# seq_file2 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_2.fasta"
# untrimmed_file = OUT_ORTSEQ_TMP / f"Untrimmed_{SAMPLE_ID}.fasta"
# trim1 = OUT_ORTSEQ_TMP / f"trimmed_orientated_all_{SAMPLE_ID}.fasta"
# trim2 = OUT_ORTSEQ_TMP / f"Filtered_J_{SAMPLE_ID}.fasta"
# trim3 = OUT_ORTSEQ_TMP / f"Filtered_reduced_{SAMPLE_ID}.fasta"
# Fail_file = OUT_FASTQ / f"Fail_filtered_{SAMPLE_ID}.fasta"
# primer_tag_file = (
#     OUT_ORTSEQ_TMP / f"Barcode_filtering_information_{SAMPLE_ID}.txt"
# )
# primer_tag_file_count = OUT_ORTSEQ_TMP / f"All_barcodes_{SAMPLE_ID}.txt"
# filtered_out = OUT_ORTSEQ / f"Filtered_ORFs_sequences_all_{SAMPLE_ID}.fasta"
# nn_orf_filtered = (
#     OUT_ORTSEQ / f"Nucleotide_ORF_filtered_all_{SAMPLE_ID}.fasta"
# )
tmp_file_orf = OUT_ORTSEQ_TMP / f"Blast_matching_{SAMPLE_ID}"
filtering_report = OUT_ORTSEQ / f"Filtering_report_{SAMPLE_ID}.txt"
# Files for clustering
att_file = OUT_NET / f"Vertex_relations_{SAMPLE_ID}.txt"
file_vertex = OUT_NET / f"Att_{SAMPLE_ID}.txt"
file_edges = OUT_NET / f"Edges_{SAMPLE_ID}.txt"
cluster_file = OUT_NET / f"Cluster_sample_identities_{SAMPLE_ID}.txt"
Reduced_file = OUT_NET / f"Fully_reduced_{SAMPLE_ID}.fasta"
checked_edges = OUT_NET / f"Checked_edges_{SAMPLE_ID}.txt"
plot_sample_ids_file = OUT_NET / f"Plot_sample_ids_{SAMPLE_ID}.txt"
file_seqs = OUT_NET / f"Sequences_{SAMPLE_ID}.txt"
tmp_file0 = OUT_NET / f"Decon_0_{SAMPLE_ID}.txt"
tmp_pre = OUT_NET / f"Pre_tmp_{SAMPLE_ID}"
tmp_file_1 = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.1"
# tmp_file = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}."


# Commands
if command_source.count("1") != 0:
    from isotyper.qualitycontrol import preliminary

    preliminary.main()

if command_source.count("2") != 0:
    from isotyper.qualitycontrol import qualitycontrol

    qualitycontrol.main()


# Clustering reads
if command_source.count("3") != 0:
    generate_networks(
        FILTERED_OUT_NT, tmp_file_1, EDGE_LENGTHS, tmp_file, att_file
    )
    deconvolute_edges(
        FILTERED_OUT_NT,
        att_file,
        file_vertex,
        file_seqs,
        tmp_file0,
        file_edges,
        READ_NUMBER_DIVISION,
    )
    get_network_clusters(file_vertex, file_edges, cluster_file)
    reduce_identical_sequences(Reduced_file, file_vertex, READ_NUMBER_DIVISION)
    get_network_input(cluster_file, plot_ids_file, file_edges, checked_edges)
