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
from typing import Union, Tuple

from isotyper.utilities import *


def reduce_sequences(trim2, trim3, primer_file):
    """Summary

    Parameters
    ----------
    trim2 : TYPE
        Description
    trim3 : TYPE
        Description
    primer_file : TYPE
        Description
    """
    const = "TRUE"
    minl = 185  # change for shorter runs
    if const == "TRUE":
        fh = open(trim3, "w")
        fh.close()
        fh = open(trim2, "r")
        seqs, head = Tree(), ""
        for header, seq in fasta_iterator(fh):
            if len(seq) >= minl:
                seq = seq.upper()
                seqs[seq][header].value = 1
                head = header.split("|")[1]
            else:
                print(seq)
        fh.close()
        out, ind, times = "", 0, len(head.split("_"))
        print("number of chains", times)
        for seq in seqs:
            f = [0] * times
            for id in seqs[seq]:
                f1 = list(map(int, id.split("__")[1].split("|")[0].split("_")))
                f = list(map(add, f, f1))
            header = (
                ">"
                + id.split("__")[0]
                + "__"
                + "_".join(map(str, f))
                + "|"
                + head
            )
            out = out + header + "\n" + seq + "\n"
            ind = ind + 1
            if ind > 500:
                write_out(out, trim3)
                out, ind = "", 0
        write_out(out, trim3)
    else:
        command1 = "cp {} {}".format(trim2, trim3)
        os.system(command1)


def get_match(primer, seq):
    """Summary

    Parameters
    ----------
    primer : TYPE
        Description
    seq : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    loc = []
    if seq.count(primer) != 0:
        for m in re.finditer(primer, seq):
            loc = loc + [m.start()]
            # loc = seq.index(primer)
    return loc


def filter_igj_genes(
    trim1, trim2, refj, primer_file, ref_const, primer_tag_file_count
):
    """Summary

    Parameters
    ----------
    trim1 : TYPE
        Description
    trim2 : TYPE
        Description
    refj : TYPE
        Description
    primer_file : TYPE
        Description
    ref_const : TYPE
        Description
    primer_tag_file_count : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    mode = "WITHIN"
    fh = open(trim2, "w")
    fh.close()
    out, ind, batch, batch_size = "", 0, 0, 100
    seqs, indent = {}, 120
    fh = open(trim1, "r")
    e_value = 10
    c = 0
    for header, sequence in fasta_iterator(fh):
        # inf = len(sequence) - indent
        inf = 0
        # if inf < 0:
        #     10
        out = out + ">" + header + "\n" + sequence[inf : len(sequence)] + "\n"
        seqs[header] = sequence
        ind, batch = ind + 1, batch + 1
        c = c + 1
        if batch >= batch_size:
            blast_match_j(out, seqs, trim1, trim2, refj, e_value)
            out, batch, seqs = "", 0, {}
    fh.close()
    if len(out) > 2:
        blast_match_j(out, seqs, trim1, trim2, refj, e_value)
        out, batch = "", 0


def check_fasta_not_empty(fh):
    """Summary

    Parameters
    ----------
    fh : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    pfh = 0
    for l in fh:
        print(l)
        if len(l) != 0:
            pfh = 1
        break
    return pfh


def get_consensus_sequence(u_seq, u_freq, tmp_file):
    """Summary

    Parameters
    ----------
    u_seq : TYPE
        Description
    u_freq : TYPE
        Description
    tmp_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    out = ""
    for i in range(0, len(u_seq)):
        out = out + ">" + str(i) + "\n" + u_seq[i] + "\n"
    fh = open(tmp_file + "txt", "w")
    fh.write(out)
    fh.close()
    if len(u_seq) > 2000:
        command1 = (
            "mafft --retree 2 --parttree "
            + "{}txt".format(tmp_file)
            + " > {}aligned".format(tmp_file)
        )
    else:
        command1 = (
            "mafft --retree 2 "
            + "{}txt".format(tmp_file)
            + " > {}aligned".format(tmp_file)
        )
    os.system(command1)
    fh = open(tmp_file + "aligned", "r")
    max_seqs = {}
    pfh = check_fasta_not_empty(fh)
    fh.close()
    if pfh == 1:
        fh = open(tmp_file + "aligned", "r")
        for header, sequence in fasta_iterator(fh):
            max_seqs[sequence.upper()] = int(header)
        fh.close()
        bases = ["A", "T", "G", "C", "-"]
        base_dict = {}
        for b in range(0, len(bases)):
            base_dict[bases[b]] = b
        consensus = ""
        start = 0
        for i in range(0, len(sequence)):
            f = [0] * len(bases)
            for s in max_seqs:
                f[base_dict[s[i]]] = f[base_dict[s[i]]] + u_freq[max_seqs[s]]
            if f[4] == 0:
                start == 1
            if max(f) * 1.0 / sum(f) >= THRESHOLD_BARCODE:
                if bases[f.index(max(f))] != "-":
                    consensus = consensus + bases[f.index(max(f))]
            else:
                if start == 1:
                    for j in range(0, 5):
                        if f[j] != 0:
                            consensus = (
                                consensus
                                + "|"
                                + bases[j]
                                + ":"
                                + str(
                                    "%s" % float("%.3g" % (f[j] * 1.0 / sum(f)))
                                )
                            )
                    consensus = consensus + "_"
                else:
                    f = f[0:4]
                    if bases[f.index(max(f))] != "-":
                        consensus = consensus + bases[f.index(max(f))]
                    else:
                        for j in range(0, 4):
                            if f[j] != 0:
                                consensus = (
                                    consensus
                                    + "|"
                                    + bases[j]
                                    + ":"
                                    + str(
                                        "%s"
                                        % float("%.3g" % (f[j] * 1.0 / sum(f)))
                                    )
                                )
                            consensus = consensus + "_"
        # else:print "ERROR", u_freq, u_seq
    else:
        consensus = ""
    return consensus


def trim_sequences_bcr_tcr(
    tmp_Tmp_file,
    Fail_file,
    Output_trim,
    gene,
    paired,
    species,
    primer_file,
    primer_tag_file,
    tmp_file,
    primer_tag_file_count,
    sample,
    ref_const,
    reverse_primer_group,
):
    """Summary

    Parameters
    ----------
    tmp_Tmp_file : TYPE
        Description
    Fail_file : TYPE
        Description
    Output_trim : TYPE
        Description
    gene : TYPE
        Description
    paired : TYPE
        Description
    species : TYPE
        Description
    primer_file : TYPE
        Description
    primer_tag_file : TYPE
        Description
    tmp_file : TYPE
        Description
    primer_tag_file_count : TYPE
        Description
    sample : TYPE
        Description
    ref_const : TYPE
        Description
    reverse_primer_group : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    forward, reverse, barcoded_j, barcoded_v, v_ref = get_primers_split(
        primer_file
    )
    fh_out = open(Output_trim, "w")
    fh_out.close()
    fh_out = open(Fail_file, "w")
    fh_out.close()
    if barcoded_j == 1 and barcoded_v == 0:
        print("J barcoded")
        single_j_barcoded_trimming_clustered(
            forward,
            reverse,
            barcoded_j,
            barcoded_v,
            tmp_Tmp_file,
            Fail_file,
            Output_trim,
            primer_tag_file,
            tmp_file,
            gene,
            paired,
            species,
            primer_file,
            primer_tag_file_count,
            ref_const,
            v_ref,
        )
    return ()


def get_primers_split(primer_file: Path) -> Tuple:
    """Summary

    Parameters
    ----------
    primer_file : Path
        path to primer file.

    Returns
    -------
    Tuple
        output forward, reverse, barcoded_j, barcoded_v and v_ref
    """
    fh = open(primer_file, "r")
    forward, reverse, v_ref = [], [], []
    barcoded_j, barcoded_v = 0, 0
    word_size = 8
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            if len(l) > 2:
                header, sequence, uni, BC, gene_specific = (
                    l[0],
                    l[1],
                    l[2],
                    l[3],
                    l[4],
                )
                gene_specific = gene_specific.upper()
                words = []
                for i in range(word_size, len(gene_specific)):
                    words.append(
                        [gene_specific[i - word_size : i], i - word_size]
                    )
                if (
                    header.count("J") != 0
                    or header.count("REV_CONST") != 0
                    or header.count("REVERSE") != 0
                ):
                    if sequence.count("N") != 0:
                        barcoded_j = 1
                    sequence = sequence.upper()
                    if header.count("HOUSEKEEPING") != 0:
                        clas = "HOUSEKEEPING"
                    else:
                        clas = "IMMUME_REC"
                    l = len(gene_specific)
                    reverse = reverse + [
                        [sequence, clas, uni, BC, gene_specific, header, words]
                    ]
                else:
                    if sequence.count("N") != 0:
                        barcoded_v = 1
                    sequence = sequence.upper()
                    if header.count("HOUSEKEEPING") != 0:
                        clas = "HOUSEKEEPING"
                    else:
                        clas = "IMMUNE_REC"
                    l = len(gene_specific)
                    forward = forward + [
                        [sequence, clas, uni, BC, gene_specific, header, words]
                    ]
            else:
                header, sequence = l[0], l[1]
                sequence = sequence.upper()
                words = []
                for i in range(word_size, len(sequence)):
                    words.append([sequence[i - word_size : i], i - word_size])
                if (
                    header.count("J") != 0
                    or header.count("REV_CONST") != 0
                    or header.count("REVERSE") != 0
                ):
                    if sequence.count("N") != 0:
                        barcoded_j = 1
                    sequence = sequence.upper()
                    clas = "IMMUME_REC"
                    reverse = reverse + [
                        [sequence, clas, sequence, header, words]
                    ]
                elif header.count("REF") == 0:
                    if sequence.count("N") != 0:
                        barcoded_v = 1
                    sequence = sequence.upper()
                    clas = "IMMUNE_REC"
                    forward = forward + [[sequence, clas, header, words]]
                else:
                    sequence = sequence.upper()
                    clas = "IMMUNE_REC"
                    v_ref = v_ref + [[sequence, clas, header, words]]
    fh.close()
    return (forward, reverse, barcoded_j, barcoded_v, v_ref)


def separate_sequences(
    primer_tag_file_count, primer_tag_file, Output_trim, ref_const
):
    """Summary

    Parameters
    ----------
    primer_tag_file_count : TYPE
        Description
    primer_tag_file : TYPE
        Description
    Output_trim : TYPE
        Description
    ref_const : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(ref_const, "r")
    word_length = 8
    word_dict = {}
    for header, sequence in fasta_iterator(fh):
        header, sequence = (
            header.replace("|", ",").split("*")[0],
            sequence.upper(),
        )
        words = []
        for i in range(word_length, len(sequence)):
            words.append([sequence[i - word_length : i], i - word_length])
        word_dict[header] = words
    fh.close()
    region_primer_type = {}
    fh = open(primer_tag_file_count, "r")
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            const_reg, seq_type, rev_primer = l[len(l) - 1], l[4], l[6]
            region_primer_type[l[0].split(":")[0]] = [
                const_reg,
                seq_type,
                rev_primer,
            ]
    fh.close()
    fh = open(primer_tag_file, "r")
    seqs = Tree()
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            if l[6] == "YES":
                consensus = l[7]
                id = l[0].split("__")[0]
                seqs[consensus][id].value = 1
    fh.close()
    print(len(seqs))
    out, ind = "", 0
    seq_info, seq_uniq = {}, Tree()
    classes = {}
    for s in seqs:
        for id in seqs[s]:
            break
        scores = assess_gene_score(s, word_dict)
        if len(scores) > 0:
            if (
                scores[0][0] > threshold
            ):  # for non-housekeeping genes, trimmed sequences to remove constant region
                if len(scores) == 1:
                    if len(scores[0][2]) > 50:
                        seq = scores[0][2]
                    else:
                        seq = s
                    seq_uniq[s][id].value = 1
                    seq_info[id] = [len(seqs[s]), seq, scores[0][1]]
                    classes[scores[0][1]] = 1
    del seqs
    print(len(classes))
    classes_all = []
    for c in classes:
        classes_all.append(c)
    classes_all.sort()
    classes_order = {}
    for i in range(0, len(classes_all)):
        classes_order[classes_all[i]] = i
    header = "_".join(classes_all)
    for s in seq_uniq:
        if len(seq_uniq[s]) > 0:
            f = [0] * len(classes_all)
            for id in seq_uniq[s]:
                c = seq_info[id][2]
                f[classes_order[c]] = f[classes_order[c]] + seq_info[id][0]
            out = (
                out
                + ">"
                + id.split("__")[0]
                + "__"
                + "_".join(map(str, f))
                + "|"
                + header
                + "\n"
                + s
                + "\n"
            )
            ind = ind + 1
            if ind > 100:
                write_out(out, Output_trim)
                out, ind = "", 0
    write_out(out, Output_trim)
    return ()


def assess_gene_score(consensus, word_dict):
    """Summary

    Parameters
    ----------
    consensus : TYPE
        Description
    word_dict : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    scores, max_score = [], 0
    for hk in word_dict:
        w = word_dict[hk]
        score = [i for i in range(len(w)) if consensus.count(w[i][0]) != 0]
        if len(score) > THRESHOLD_GENE_SCORE:
            if len(score) >= max_score:
                consensus1 = consensus
                if (
                    hk.count("HOUSEKEEPING") == 0
                    and len(score) > THRESHOLD_GENE_SCORE
                ):
                    start = [
                        consensus.index(w[score[i]][0]) - w[score[i]][1]
                        for i in range(len(score))
                    ]
                    start = max(set(start), key=start.count)
                    consensus1 = consensus[0:start]
                if max_score == len(score):
                    scores.append([len(score), hk, consensus1])
                else:
                    scores = []
                    scores.append([len(score), hk, consensus1])
                max_score = len(score)
    return scores


def single_j_barcoded_trimming_clustered(
    forward,
    reverse,
    barcoded_j,
    barcoded_v,
    tmp_Tmp_file,
    Fail_file,
    Output_trim,
    primer_tag_file,
    tmp_file,
    gene,
    paired,
    species,
    primer_file,
    primer_tag_file_count,
    ref_const,
    v_ref,
):
    """Summary

    Parameters
    ----------
    forward : TYPE
        Description
    reverse : TYPE
        Description
    barcoded_j : TYPE
        Description
    barcoded_v : TYPE
        Description
    tmp_Tmp_file : TYPE
        Description
    Fail_file : TYPE
        Description
    Output_trim : TYPE
        Description
    primer_tag_file : TYPE
        Description
    tmp_file : TYPE
        Description
    gene : TYPE
        Description
    paired : TYPE
        Description
    species : TYPE
        Description
    primer_file : TYPE
        Description
    primer_tag_file_count : TYPE
        Description
    ref_const : TYPE
        Description
    v_ref : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    read_untrimmed_file_single(
        tmp_Tmp_file,
        Fail_file,
        Output_trim,
        gene,
        paired,
        species,
        primer_file,
        primer_tag_file,
        tmp_file,
        primer_tag_file_count,
        forward,
        reverse,
        v_ref,
    )
    check_barcodes_malbac(
        primer_tag_file_count,
        primer_tag_file,
        Fail_file,
        Output_trim,
    )
    separate_sequences(
        primer_tag_file_count, primer_tag_file, Output_trim, ref_const
    )
    return ()


def check_barcodes_malbac(
    primer_tag_file_count, primer_tag_file, Fail_file, Output_trim
):
    """Summary

    Parameters
    ----------
    primer_tag_file_count : TYPE
        Description
    primer_tag_file : TYPE
        Description
    Fail_file : TYPE
        Description
    Output_trim : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    header_txt = (
        "#ID\tnumber_of_reads\ttotal_reads_with_BC\tJ_barcode\tV_barcode\t"
        + "bp_mismatches_to_consensus\tBC_accepted\tconsensus\tsequence\n"
    )
    outs = ["", header_txt]
    files = [Output_trim, primer_tag_file]
    for i in range(0, len(files)):
        fh = open(files[i], "w")
        fh.write(outs[i])
        fh.close()
    fh = open(primer_tag_file_count, "r")
    seqs = Tree()
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            v_tag, j_tag, sequence, header = l[2], l[1], l[3], l[0]
            seqs[j_tag + "\t" + v_tag][sequence][header].value = 1
    fh.close()
    print(len(seqs), "Unique tags")
    total_tags, ind, passed_seqs_total = 0, 0, 0
    fail_less_than_threshold = 0
    outp = ""
    for t1 in seqs:
        u_seq, u_freq, u_header = [], [], []
        for s in seqs[t1]:
            f = 0
            for h in seqs[t1][s]:
                f = f + int(h.split(":")[1])
            total_tags, u_seq, u_freq, u_header, ind = (
                total_tags + f,
                u_seq + [s],
                u_freq + [f],
                u_header + [h],
                ind + 1,
            )
        f = sum(u_freq)
        h = h.split(":")[0].split("__")[0] + "__" + str(sum(u_freq))
        if sum(u_freq) > 20:
            print("BC group size:\t", len(u_freq), t1.split())
        if len(u_freq) == 1:
            passed_seqs_total = passed_seqs_total + f
            outp = (
                outp
                + h
                + "\t"
                + str(f)
                + "\t"
                + str(f)
                + "\t"
                + t1
                + "\t0\tYES\t"
                + s
                + "\t"
                + s
                + "\n"
            )
        elif len(u_freq) < 500:
            if max(u_freq) > sum(u_freq) * THRESHOLD_BARCODE:
                nz = [i for i in range(len(u_freq)) if u_freq[i] != max(u_freq)]
                passed_seqs_total = passed_seqs_total + f
                consensus = u_seq[nz[0]]
                outp = (
                    outp
                    + "\t".join(
                        map(
                            str,
                            [
                                h,
                                len(u_freq),
                                sum(u_freq),
                                t1,
                                0,
                                "YES",
                                consensus,
                                consensus,
                            ],
                        )
                    )
                    + "\n"
                )
            elif len(u_freq) > 15:  # clustering first then alignment
                out_cluster, ids = "", {}
                for i in range(0, len(u_seq)):
                    out_cluster = (
                        out_cluster + ">" + u_header[i] + "\n" + u_seq[i] + "\n"
                    )
                    ids[u_header[i]] = [u_seq[i], u_freq[i]]
                consensus, pass_consensus = get_consensus_sequence_large(
                    out_cluster,
                    Fail_file,
                    len(u_seq[i]),
                    ids,
                    sum(u_freq),
                    tmp_file,
                )
                if consensus != "" and pass_consensus == 1:
                    outp = (
                        outp
                        + "\t".join(
                            map(
                                str,
                                [
                                    h,
                                    len(u_freq),
                                    sum(u_freq),
                                    t1,
                                    0,
                                    "YES",
                                    consensus,
                                    consensus,
                                ],
                            )
                        )
                        + "\n"
                    )
                    passed_seqs_total = passed_seqs_total + f
                else:
                    fail_less_than_threshold = fail_less_than_threshold + 1
            else:
                consensus, pass_consensus = get_consensus_sequence_cluster(
                    u_seq, u_freq, tmp_file
                )
                if consensus.count("_") == 0 and pass_consensus == 1:
                    outp = (
                        outp
                        + "\t".join(
                            map(
                                str,
                                [
                                    h,
                                    len(u_freq),
                                    sum(u_freq),
                                    t1,
                                    0,
                                    "YES",
                                    consensus,
                                    consensus,
                                ],
                            )
                        )
                        + "\n"
                    )
                    passed_seqs_total = passed_seqs_total + f
                else:
                    fail_less_than_threshold = fail_less_than_threshold + 1
        if ind > 200:
            write_out(outp, primer_tag_file)
            outp, ind = "", 0
    write_out(outp, primer_tag_file)
    outp, ind = "", 0
    print(total_tags, passed_seqs_total, fail_less_than_threshold)
    return ()


def get_consensus_sequence_cluster(u_seq, u_freq, tmp_file):
    """Summary

    Parameters
    ----------
    u_seq : TYPE
        Description
    u_freq : TYPE
        Description
    tmp_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    out = ""
    for i in range(0, len(u_seq)):
        out = out + ">" + str(i) + "\n" + u_seq[i] + "\n"
    fh = open(tmp_file + "txt", "w")
    fh.write(out)
    fh.close()
    if len(u_seq) > 2000:
        # command1 = "mafft --retree 2 " + insert + " " + tmp_file + "txt > " + tmp_file + "aligned"
        command1 = (
            "mafft --retree 2 --parttree "
            + "{}txt".format(tmp_file)
            + " > {}aligned".format(tmp_file)
        )
    else:
        command1 = (
            "mafft --retree 2 "
            + "{}txt".format(tmp_file)
            + " > {}aligned".format(tmp_file)
        )
    os.system(command1)
    fh = open(tmp_file + "aligned", "r")
    max_seqs = {}
    for header, sequence in fasta_iterator(fh):
        max_seqs[sequence.upper()] = int(header)
    fh.close()
    bases = ["A", "T", "G", "C", "-"]
    base_dict = {}
    for b in range(0, len(bases)):
        base_dict[bases[b]] = b
    consensus = ""
    start, pass_consensus = 0, 1
    for i in range(0, len(sequence)):
        f = [0] * len(bases)
        for s in max_seqs:
            if s[i] not in base_dict:
                print(i, s[i], max_seqs)
            f[base_dict[s[i]]] = f[base_dict[s[i]]] + u_freq[max_seqs[s]]
        if f[4] == 0:
            start == 1
        if max(f) * 1.0 / sum(f) >= THRESHOLD_BARCODE:
            if bases[f.index(max(f))] != "-":
                consensus = consensus + bases[f.index(max(f))]
        else:
            pass_consensus = 0
            if start == 1:
                for j in range(0, 5):
                    if f[j] != 0:
                        consensus = (
                            consensus
                            + "|"
                            + bases[j]
                            + ":"
                            + str("%s" % float("%.3g" % (f[j] * 1.0 / sum(f))))
                        )
                consensus = consensus + "_"
            else:
                f = f[0:4]
                if bases[f.index(max(f))] != "-":
                    consensus = consensus + bases[f.index(max(f))]
                else:
                    for j in range(0, 4):
                        if f[j] != 0:
                            consensus = (
                                consensus
                                + "|"
                                + bases[j]
                                + ":"
                                + str(
                                    "%s" % float("%.3g" % (f[j] * 1.0 / sum(f)))
                                )
                            )
                        consensus = consensus + "_"
    return (consensus, pass_consensus)


def get_consensus_sequence_large(
    out_cluster, Fail_file, l1, ids, sum_u_freq, tmp_file
):
    """Summary

    Parameters
    ----------
    out_cluster : TYPE
        Description
    Fail_file : TYPE
        Description
    l1 : TYPE
        Description
    ids : TYPE
        Description
    sum_u_freq : TYPE
        Description
    tmp_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(Fail_file, "w")
    fh.write(out_cluster)
    fh.close()
    cluster_i(Fail_file, Fail_file + "cls", (l1 - 3.0) / l1)
    cluster = Tree()
    fh = open(Fail_file + "cls.bak.clstr", "r")
    for l in fh:
        l = l.strip().split()
        clust, id = l[0], l[2].replace("...", "").replace(">", "")
        cluster[clust][id].value = 1
    fh.close()
    max_s, max_clust = 0, ""
    for c in cluster:
        f = 0
        for id in cluster[c]:
            f = f + ids[id][1]
        if max_s < f:
            max_clust, max_s = c, f
    consensus, pass_consensus = "", 0
    if sum_u_freq * THRESHOLD_BARCODE < max_s:
        seq_align, s_freq = [], []
        for id in cluster[max_clust]:
            seq_align.append(ids[id][0])
            s_freq.append(ids[id][1])
        consensus = get_consensus_sequence(seq_align, s_freq, tmp_file)
        if consensus.count("_") == 0 and len(consensus) > 3:
            pass_consensus = 1
        else:
            pass_consensus = 0
    return (consensus, pass_consensus)


def read_untrimmed_file_single(
    tmp_Tmp_file,
    Fail_file,
    Output_trim,
    gene,
    paired,
    species,
    primer_file,
    primer_tag_file,
    tmp_file,
    primer_tag_file_count,
    forward,
    reverse,
    v_ref,
):
    """Summary

    Parameters
    ----------
    tmp_Tmp_file : TYPE
        Description
    Fail_file : TYPE
        Description
    Output_trim : TYPE
        Description
    gene : TYPE
        Description
    paired : TYPE
        Description
    species : TYPE
        Description
    primer_file : TYPE
        Description
    primer_tag_file : TYPE
        Description
    tmp_file : TYPE
        Description
    primer_tag_file_count : TYPE
        Description
    forward : TYPE
        Description
    reverse : TYPE
        Description
    v_ref : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    for f in [primer_tag_file_count]:
        fh = open(f, "w")
        fh.close()
    fh = open(tmp_Tmp_file, "r")
    minl, maxl = 110, 1000000
    if gene == "HEAVY" or gene == "IGH":
        minl = 120
    if gene == "KAPPA" or gene == "IGK":
        minl = 110
    if gene == "LAMBDA" or gene == "IGL":
        (minl, maxl) = (90, 150)
    seqs1, t = Tree(), 0
    for header, seq in fasta_iterator(fh):
        seq = seq.upper()
        seqs1[seq][header].value = 1
        t = t + 1
    out, ind = (
        "#ID\tJ_tag\tV_tag\tSequence\ttype_rev\ttype_for\tprimer_rev\tprimer_for\n",
        0,
    )
    total, pass_r, pass_f, pass_all = 0, 0, 0, 0
    for seq in seqs1:
        for header in seqs1[seq]:
            break
        header = header + ":" + str(len(seqs1[seq]))
        passes, j_tag, v_tag = 0, "", ""
        type_rev, type_for = "", ""
        primer_rev, primer_for = "", ""
        total = total + 1
        for i in range(0, len(reverse)):
            pj = get_match(reverse[i][2], seq)
            if max(pj + [-1]) == -1 or len(pj) == 0:
                seq = reverse_comp(seq)
                pj = get_match(reverse[i][2], seq)
            if max(pj + [-1]) != -1:
                bc_len = len(reverse[i][3])
                if pj[0] > bc_len - 3:
                    j_tag = seq[pj[0] - bc_len : pj[0]]
                    if len(j_tag) > bc_len / 2:
                        seq = seq[pj[0] + bc_len : len(seq)]
                        type_rev, primer_rev = reverse[i][1], reverse[i][5]
                        passes = 1
                        break
        if passes != 1:
            for i in range(0, len(reverse)):
                words = reverse[i][6]
                p = []
                for w in words:
                    pj = get_match(w[0], seq)
                    if max(pj + [-1]) != -1 and len(pj) > 0:
                        if pj[0] < len(seq) / 2:
                            p = p + [max([0, pj[0] - w[1]])]
                        else:
                            pj[0] = -1
                if len(p) > 1:
                    pj = max(set(p), key=p.count)
                    bc_len = len(reverse[i][3])
                    if pj > min([bc_len - 3]):
                        j_tag = seq[pj - bc_len : pj + 1]
                        if len(j_tag) > bc_len / 2:
                            seq = seq[pj + bc_len : len(seq)]
                            seq = reverse_comp(seq)
                            type_rev, primer_rev = reverse[i][1], reverse[i][5]
                            passes = 1
                            break
            if passes != 1:
                seq = reverse_comp(seq)
                for i in range(0, len(reverse)):
                    words = reverse[i][6]
                    p = []
                    for w in words:
                        pj = get_match(w[0], seq)
                        if max(pj + [-1]) != -1 and len(pj) > 0:
                            if pj[0] < len(seq) / 2:
                                p = p + [max([0, pj[0] - w[1]])]
                            else:
                                pj[0] = -1
                    if len(p) > 1:
                        bc_len = len(reverse[i][3])
                        pj = max(set(p), key=p.count)
                        if pj > min([bc_len - 3]):
                            j_tag = seq[pj - bc_len : pj + 1]
                            if len(j_tag) > bc_len / 2:
                                seq = seq[pj + bc_len : len(seq)]
                                seq = reverse_comp(seq)
                                type_rev, primer_rev = (
                                    reverse[i][1],
                                    reverse[i][5],
                                )
                                passes = 1
                                break
        if passes == 1:
            pass_r = pass_r + 1
            for i in range(0, len(forward)):
                pv = get_match(forward[i][0], seq)
                if max(pv + [-1]) != -1:
                    v_tag = "-"
                    seq = seq[pv[0] + len(forward[i][0]) - 1 : len(seq)]
                    type_for, primer_for = forward[i][1], forward[i][2]
                    passes = 2
                    break
            if passes != 2:
                for i in range(0, len(forward)):
                    words = forward[i][3]
                    p = 0
                    for w in words:
                        pv = get_match(w[0], seq)
                        if max(pv + [-1]) != -1 and len(pv) > 0:
                            if pv[0] < len(seq) / 2:
                                p = p + 1
                            else:
                                pv[0] = -1
                    if max(pv + [-1]) != -1 and len(pv) > 0 and p > 1:
                        v_tag = "-"
                        seq = seq[pv[0] + len(forward[i][2]) - 1 : len(seq)]
                        type_for, primer_for = forward[i][1], forward[i][2]
                        passes = 2
                        break
            if passes != 2:
                p = []
                for i in range(0, len(v_ref)):
                    words = v_ref[i][3]
                    p1 = []
                    for w in words:
                        pv = get_match(w[0], seq)
                        if max(pv + [-1]) != -1 and len(pv) > 0:
                            if pv[0] < len(seq) / 2:
                                p = p + [pv[0]]
                    if len(p1) > 0:
                        p1 = p1 + [min(p1)]
                if len(p) >= 2:
                    pv = min(p)
                    v_tag = "-"
                    seq = seq[pv : len(seq)]
                    type_for, primer_for = v_ref[i][1], v_ref[i][2]
                    passes = 2
        if passes == 2:
            pass_f = pass_f + 1
            if len(seq) > minl and len(seq) < maxl:
                if seq.count("N") == 0:
                    # if(type_rev==type_for):
                    out = (
                        out
                        + "\t".join(
                            map(
                                str,
                                [
                                    header,
                                    j_tag,
                                    v_tag,
                                    seq,
                                    type_rev,
                                    type_for,
                                    primer_rev,
                                    primer_for,
                                ],
                            )
                        )
                        + "\n"
                    )
                    pass_all = pass_all + 1
                    ind = ind + 1
                    if ind > 200:
                        write_out(out, primer_tag_file_count)
                        out, ind = "", 0
    fh.close()
    write_out(out, primer_tag_file_count)
    out, ind = "", 0
    print(
        "\tTotal sequences:\t"
        + str(total)
        + "\tNumber of REV primers found:\t"
        + str(pass_r)
        + "\tNumber of FOR primers found:\t"
        + str(pass_f)
        + "\tPassed sequences:\t"
        + str(pass_all)
    )
    return ()


def get_sequences(file):
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
    command = "gunzip {}".format(file)
    os.system(command)
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        header = header.replace(":", "").split()[0].replace("-", "")
        header = header.split("#")[0].split("/")[0]
        seqs[header] = sequence
        # if(len(seqs)>100000):break ####### remove
    fh.close()
    return seqs


def trim(s1, s2, l1, l2, indent, length):
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
    indent : TYPE
        Description
    length : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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


def join_reads(s1, s2, length):
    """Summary

    Parameters
    ----------
    s1 : TYPE
        Description
    s2 : TYPE
        Description
    length : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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


def get_paired_reads_overlapping(
    file1, file2, outfile, gene, paired, id, method
):
    """Summary

    Parameters
    ----------
    file1 : TYPE
        Description
    file2 : TYPE
        Description
    outfile : TYPE
        Description
    gene : TYPE
        Description
    paired : TYPE
        Description
    id : TYPE
        Description
    method : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    seqs1 = get_sequences(file1)
    seqs2 = get_sequences(file2)
    print("Forward reads:", len(seqs1), "Reverse reads:", len(seqs2))
    out, ind = "", 0
    fh = open(outfile, "w")
    fh.close()
    tot, f_joining, total, length, fail_no_pair = 0, 0, 0, 30, 0
    if method in ["5RACE"]:
        length = 20
    if gene == "LAMBDA" or gene == "KAPPA" or gene == "IGK":
        length = 30
    for id in seqs1:
        if id in seqs2:
            total = total + 1
            (seq, failed) = join_reads(seqs1[id], seqs2[id], length)
            if failed == 0 and len(seq) > 180:
                ind, tot = ind + 1, tot + 1
                out = out + ">" + id + "\n" + seq + "\n"
            else:
                seq2 = reverse_comp(seqs2[id])
                (seq, failed) = join_reads(seqs1[id], seq2, length)
                if failed == 0 and len(seq) > 180:
                    ind, tot = ind + 1, tot + 1
                    out = out + ">" + id + "\n" + seq + "\n"
                else:
                    f_joining = f_joining + 1
                    # print seqs1[id], seqs2[id]
            if ind > 100:
                write_out(out, outfile)
                out, ind = "", 0
        else:
            fail_no_pair = fail_no_pair + 1
    write_out(out, outfile)
    del out, ind
    print(
        id + "\tTotal sequences:",
        total,
        "\tNumber of Failed sequences:",
        f_joining,
        "\tPercentage sequences failed:",
        f_joining * 100.0 / total,
        "%\tTotal remaining:",
        tot,
        "\tNo pairs:",
        fail_no_pair,
    )
    return ()


def calculate_orf_length(codon, sequence, type, gene):
    """Summary

    Parameters
    ----------
    codon : TYPE
        Description
    sequence : TYPE
        Description
    type : TYPE
        Description
    gene : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    three_frames = [sequence[i:] for i in range(3)]
    three_frames_translated = [
        translate(frame, codon) for frame in three_frames
    ]
    accepted = []
    # found = 0
    for i in range(0, len(three_frames_translated)):
        frame = three_frames_translated[i]
        if frame.count("-") == 0:
            accepted.append((frame, i))
    accepts = 1
    if len(accepted) == 0:
        accepts = 0
    return (accepted, accepts)


def get_sequences_ref(file, word):
    """Summary

    Parameters
    ----------
    file : TYPE
        Description
    word : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(file, "r")
    dict = Tree()
    for header, sequence in fasta_iterator(fh):
        for i in range(1, len(sequence) - word - 1):
            dict[sequence[i : (i + word)]][header].value = 1
    fh.close()
    return dict


def get_max_orf(ORF, word, dict):
    """Summary

    Parameters
    ----------
    ORF : TYPE
        Description
    word : TYPE
        Description
    dict : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    score = []
    codon = []
    if len(ORF) > 1:
        for si in ORF:
            s = si[0]
            sc = 0
            for i in range(1, len(s) - word):
                if s[i : (i + word)] in dict:
                    sc = sc + 1
            score.append(sc)
            codon.append(si[1])
        maxi = score.index(max(score))
        orf = ORF[maxi][0]
        mscore = max(score)
        codon = ORF[maxi][1]
    else:
        for si in ORF:
            orf = si[0]
            mscore = 1000
            codon = si[1]
    return (orf, mscore, codon)


def blast_match_j(out, seqs, trim1, trim2, refj, e_value):
    """Summary

    Parameters
    ----------
    out : TYPE
        Description
    seqs : TYPE
        Description
    trim1 : TYPE
        Description
    trim2 : TYPE
        Description
    refj : TYPE
        Description
    e_value : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(trim1 + "_blast_J", "w")
    fh.write(out)
    fh.close()
    command1 = (
        bin_path
        + "blastall -p blastn -a 10 -d {}".format(refj)
        + " -e {}".format(str(e_value))
        + " -i {}_blast_J".format(trim1)
        + " -o {}_blast_J_results".format(trim1)
        + " -b 1 -m 8 -W 4"
    )
    os.system(command1)
    fh = open(trim1 + "_blast_J_results", "r")
    out, done = "", {}
    for l in fh:
        l = l.strip().split()
        if l[0] not in done:
            out = out + ">" + l[0] + "\n" + seqs[l[0]] + "\n"
            done[l[0]] = 1
    fh.close()
    fh = open(trim2, "a")
    fh.write(out)
    fh.close()
    return ()


def orf_calculation_single(
    Output_trim,
    Filtered_out1,
    nn_orf_filtered,
    dir_ind,
    gene,
    ref,
    refj,
    ref_protein,
    refjp,
    tmp_file,
):
    """Summary

    Parameters
    ----------
    Output_trim : TYPE
        Description
    Filtered_out1 : TYPE
        Description
    nn_orf_filtered : TYPE
        Description
    dir_ind : TYPE
        Description
    gene : TYPE
        Description
    ref : TYPE
        Description
    refj : TYPE
        Description
    ref_protein : TYPE
        Description
    refjp : TYPE
        Description
    tmp_file : TYPE
        Description
    """
    get_protein_sequences(
        Output_trim,
        Filtered_out1,
        nn_orf_filtered,
        dir_ind,
        gene,
        ref,
        refj,
        ref_protein,
        refjp,
        tmp_file,
    )
    get_nucleotide_sequences(
        Output_trim,
        Filtered_out1,
        nn_orf_filtered,
        dir_ind,
        gene,
        ref,
        refj,
        ref_protein,
        refjp,
        tmp_file,
    )


def get_protein_sequences(
    Output_trim,
    Filtered_out1,
    nn_orf_filtered,
    dir_ind,
    gene,
    ref,
    refj,
    ref_protein,
    refjp,
    tmp_file,
):
    """Summary

    Parameters
    ----------
    Output_trim : TYPE
        Description
    Filtered_out1 : TYPE
        Description
    nn_orf_filtered : TYPE
        Description
    dir_ind : TYPE
        Description
    gene : TYPE
        Description
    ref : TYPE
        Description
    refj : TYPE
        Description
    ref_protein : TYPE
        Description
    refjp : TYPE
        Description
    tmp_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(Filtered_out1, "w")
    fh.close()
    fh = open(nn_orf_filtered, "w")
    fh.close()
    (codon) = get_codons()
    # (word, v_match) = (4, 45)
    (word) = 4
    (dict1) = get_sequences_ref(ref_protein, word)
    fh = open(Output_trim, "r")
    number_seqs, number_accepted, batch, indb = 0, 0, 1000, 0
    # V_region, seqs, ORFa, ORFb = '', Tree(), {}, {}
    seqs, ORFa = Tree(), {}
    batch_number, orf_fail, p_seq = 0, 0, {}
    out = ""
    for header, seq in fasta_iterator(fh):
        if seq.count("N") == 0:
            seq = seq.replace("-", "")
            header = header.split("#")[0]
            if header.count("__") == 0:
                freq = get_freq(header)
                header = header.replace(":", "") + "__" + str(freq)
            (ORF1, accept1) = calculate_orf_length(codon, seq, "V", gene)
            ORFa[header] = ORF1
            number_seqs = number_seqs + 1  # freq
            if accept1 == 0:
                orf_fail = orf_fail + 1
                # if(int(header.split("__")[1])>100):
                #  print header, ORF1, seq
            if accept1 != 0:
                if len(ORF1) > 1:
                    min_score, found = 2, 0
                    (O1, score1, codon1) = get_max_orf(ORF1, word, dict1)
                    if score1 > min_score:
                        min_score = score1
                        p_seq[header] = O1
                        seq = seq[codon1 : len(seq)]
                        found = 1
                    if found == 1:
                        number_accepted = (
                            number_accepted + 1
                        )  # int(header.split("__")[1])
                else:
                    p_seq[header] = ORF1[0][0]
                    number_accepted = (
                        number_accepted + 1
                    )  # int(header.split("__")[1])
                    seq = seq[ORF1[0][1] : len(seq)]
                seqs[seq][header].value = 1
                indb = indb + 1
                if indb >= batch:
                    batch_number = batch_number + 1
                    out = ""
                    for id in p_seq:
                        out = out + ">" + id + "\n" + p_seq[id] + "\n"
                    write_out(out, Filtered_out1)
                    p_seq = {}
                    out = ""
                    # break ############################
    fh.close()
    write_out(out, Filtered_out1)
    out, ind = "", 0
    for id in p_seq:
        out = out + ">" + id + "\n" + p_seq[id] + "\n"
        ind = ind + 1
        if ind > 500:
            write_out(out, Filtered_out1)
            out, ind = "", 0
    write_out(out, Filtered_out1)
    return ()


def get_nucleotide_sequences(
    Output_trim,
    Filtered_out1,
    nn_orf_filtered,
    dir_ind,
    gene,
    ref,
    refj,
    ref_protein,
    refjp,
    tmp_file,
):
    """Summary

    Parameters
    ----------
    Output_trim : TYPE
        Description
    Filtered_out1 : TYPE
        Description
    nn_orf_filtered : TYPE
        Description
    dir_ind : TYPE
        Description
    gene : TYPE
        Description
    ref : TYPE
        Description
    refj : TYPE
        Description
    ref_protein : TYPE
        Description
    refjp : TYPE
        Description
    tmp_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(nn_orf_filtered, "w")
    fh.close()
    fh = open(Filtered_out1, "r")
    ids = {}
    for header, sequence in fasta_iterator(fh):
        ids[header] = 1
    fh.close()
    out, ind = "", 0
    fh = open(Output_trim, "r")
    done = {}
    for header, sequence in fasta_iterator(fh):
        if header in ids:
            out = out + ">" + header + "\n" + sequence + "\n"
            done[header] = 1
            ind = ind + 1
            if ind > 500:
                write_out(out, nn_orf_filtered)
                out, ind = "", 0
    write_out(out, nn_orf_filtered)
    fh.close()
    print(len(ids), len(done))
    return ()


def get_number_sequences(file):
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
    count = 0
    process = subprocess.Popen(
        ["wc", "-l", file], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = process.communicate()
    c = out.decode("utf-8").split()
    if c[0] not in ["wc:", "ls:"]:
        fh = open(file, "r")
        for header, sequence in fasta_iterator(fh):
            count = count + 1
        fh.close()
    else:
        count = -1
    return count


def get_reduced_number_sequences_multi_constants(file):
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
    count = 0
    process = subprocess.Popen(
        ["wc", "-l", file], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = process.communicate()
    c = out.decode("utf-8").split()
    if c[0] != "ls:" and c[0] != "0":
        fh = open(file, "r")
        for header, sequence in fasta_iterator(fh):
            count = count + sum(
                map(int, header.split("__")[1].split("|")[0].split("_"))
            )
        fh.close()
    return count


def count_barcodes(primer_tag_file):
    """Summary

    Parameters
    ----------
    primer_tag_file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(primer_tag_file, "r")
    bc_uniq, bc_count, bc_failed = {}, 0, 0
    uniq_seq = 0
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            bc = l[3] + ":" + l[4]
            if l[1] != "NA":
                freq = int(l[1])  # int(l[0].split("__")[1])
                uniq_seq = uniq_seq + 1
                bc_uniq[bc] = 1
                bc_count = bc_count + freq
                if l[6] != "YES":
                    bc_failed = bc_failed + 1
    fh.close()
    n_barcodes, uniq_sequences, total_sequences_included_before_bc = (
        len(bc_uniq),
        uniq_seq,
        bc_count,
    )
    return (n_barcodes, uniq_sequences, total_sequences_included_before_bc)


def get_read_report(
    Seq_file1,
    Seq_file2,
    tmp_Tmp_file,
    trim1,
    nn_orf_filtered,
    filtering_report,
    id,
    species,
    gene,
    dir,
    primer_tag_file_count,
    primer_file,
    method,
    barcode_group,
):
    """Summary

    Parameters
    ----------
    Seq_file1 : TYPE
        Description
    Seq_file2 : TYPE
        Description
    tmp_Tmp_file : TYPE
        Description
    trim1 : TYPE
        Description
    nn_orf_filtered : TYPE
        Description
    filtering_report : TYPE
        Description
    id : TYPE
        Description
    species : TYPE
        Description
    gene : TYPE
        Description
    dir : TYPE
        Description
    primer_tag_file_count : TYPE
        Description
    primer_file : TYPE
        Description
    method : TYPE
        Description
    barcode_group : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    barcoded_j, barcoded_v = 1, 0
    print(barcoded_j, barcoded_v)
    (
        n_barcodes,
        uniq_sequences,
        total_sequences_included_before_bc,
    ) = count_barcodes(primer_tag_file)
    print(n_barcodes, uniq_sequences, total_sequences_included_before_bc)
    raw1, raw2, joined, gene_matching = (
        get_number_sequences(Seq_file1),
        get_number_sequences(Seq_file2),
        get_number_sequences(tmp_Tmp_file),
        get_reduced_number_sequences_multi_constants(trim1),
    )
    # count_bc_found, count_uniq_bcs = -1, -1
    orf = get_reduced_number_sequences_multi_constants(nn_orf_filtered)
    number_unique_seqs = get_number_sequences(nn_orf_filtered)
    out = (
        "Directory\tSample\tSpecies\tGene\t% reads retained\tN raw reads (1)\t"
        + "N raw reads (2)\tN joined reads\tN reads with BCs\tN uniq BCs\tN reads gene matched\t"
        + "N reads w/t ORF\tUnique sequences\n"
    )
    orf_perc = str(orf * 100.0 / min([raw1, raw2]))
    if min([raw1, raw2]) == -1:
        orf_perc = "NA"
    out = (
        out
        + dir
        + "\t"
        + id
        + "\t"
        + species
        + "\t"
        + gene
        + "\t"
        + orf_perc
        + "\t"
        + str(raw1)
        + "\t"
        + str(raw2)
        + "\t"
        + str(joined)
        + "\t"
        + str(n_barcodes)
        + "\t"
        + str(uniq_sequences)
        + "\t"
        + str(gene_matching)
        + "\t"
        + str(orf)
        + "\t"
        + str(number_unique_seqs)
        + "\n"
    )
    fh = open(filtering_report, "w")
    fh.write(out)
    fh.close()
    print(out)
    return ()


def cluster_i(Reduced_file, tmp_file, diff):
    """Summary

    Parameters
    ----------
    Reduced_file : TYPE
        Description
    tmp_file : TYPE
        Description
    diff : TYPE
        Description
    """
    # cd_hit_directory = "/nfs/users/nfs_k/kt16/BCRSeq/BIN/cd-hit-v4.5.7-2011-12-16/"
    # cd_hit_directory = "/lustre/scratch117/cellgen/team297/kt16/BCRSeq/BIN/cd-hit-v4.5.7-2011-12-16/"
    # cd_hit_directory = bin_path + "cd-hit-v4.5.7-2011-12-16/"
    command = (
        "cd-hit -i "
        + Reduced_file
        + " -o "
        + tmp_file
        + " -c "
        + str(diff)
        + " -g 1 -d 180 -T 10 -M 0 -AL 40 -bak 1 -p 1"
    )
    os.system(command)


def get_seqs_single(file):
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
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split("__")[0]] = sequence
    fh.close()
    return seqs


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

    Returns
    -------
    TYPE
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
    del out
    return ()


def get_diff(s1, s2, mismatch):
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


def trim_sequences(s1, s2, l1, l2):
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


def write_out(out, file):
    """Summary

    Parameters
    ----------
    out : TYPE
        Description
    file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(file, "a")
    fh.write(out)
    fh.close()
    return ()


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

    Returns
    -------
    TYPE
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
                    (id, seqs[id.split("__")[0]], len(seqs[id.split("__")[0]]))
                )
                t = t + 1
            for c1 in coclust[c]:
                for id in cluster[c1]:
                    clust_seqs.append(
                        (
                            id,
                            seqs[id.split("__")[0]],
                            len(seqs[id.split("__")[0]]),
                        )
                    )
                    t = t + 1
            clust_seqs = sorted(clust_seqs, key=itemgetter(2), reverse=True)
            if len(clust_seqs) > 1:
                if len(clust_seqs) > 500:
                    print(len(clust_seqs), total, ind)
                get_similarity_single(clust_seqs, file_out, c)
    return ()


def get_clusters(file):
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


def get_cluster_sizes_single(file):
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


def get_vaguely_similar_seqs(s1, s2, mis):
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


def count_diffs(s1, s2, mis):
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

    Returns
    -------
    TYPE
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
    return ()


def get_coclustered(file):
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

    Returns
    -------
    TYPE
        Description
    """
    fh = open(seq_file, "r")
    freq_id = {}
    for header, sequence in fasta_iterator(fh):
        freq_id[header.split("__")[0]] = header
    fh.close()
    inverse, raw = get_inverse_ids(file_seqs, file_vertex, freq_id)
    edges, edges23 = Tree(), Tree()
    fh1 = open(att_file, "r")
    for l in fh1:
        l = l.strip()
        l1 = l.split()
        if int(l1[0]) == 1 or int(l1[0]) == 2:
            id1, id2 = l1[1].split("__")[0], l1[2].split("__")[0]
            id1, id2 = freq_id[id1], freq_id[id2]
            edges[id2][id1].value = 1
    fh1.close()
    print_single_edges(file_edges, inverse, edges, tmp_file_1, raw)
    del inverse, edges, edges23
    return ()


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
            inverse[freq_id[l[0].split("__")[0]]] = freq_id[l[1].split("__")[0]]
            raw[l[1].split("__")[0]] = 1
            ind = ind + 1
    fh.close()
    fh = open(file_vertex, "r")
    for l in fh:
        l = l.strip().split()
        if l[0].split("__")[0] in raw:
            raw[l[0].split("__")[0]] = l[0]
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

    Returns
    -------
    TYPE
        Description
    """
    fh = open(seq_file, "r")
    all, seqs, freq_id = {}, {}, {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split("__")[0]] = sequence
        all[header.split("__")[0]] = sequence
        freq_id[header.split("__")[0]] = header
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
    # same, inverse, inv, out, ind, length = {}, {}, {}, '', 0, {}
    same, inverse, out, ind, length = {}, {}, "", 0, {}
    fh = open(file_seqs, "w")
    fh1.close()
    j = header
    for c in same1:
        sub_same = Tree()
        for id1 in same1[c]:
            for id2 in same1[c][id1]:
                sub_same[id1.split("__")[0]][id2.split("__")[0]].value = 1
                sub_same[id2.split("__")[0]][id1.split("__")[0]].value = 1
        (sub_same, sub_inv) = deconvolute_same_array(sub_same)
        for i in sub_same:
            s = seqs[i.split("__")[0]]
            total = 0
            mins = s
            for j in sub_same[i]:
                j1 = freq_id[j.split("__")[0]]
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
                inverse[j.split("__")[0]] = i
                out = (
                    out
                    + ">"
                    + j
                    + "||"
                    + i
                    + "\n"
                    + seqs[j.split("__")[0]]
                    + "\n"
                )
                s = seqs[j.split("__")[0]]
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
        len(freq_id[j.split("__")[0]].split(read_number_division)[1].split("|"))
        >= 2
    ):
        info = (
            "|"
            + freq_id[j.split("__")[0]]
            .split(read_number_division)[1]
            .split("|")[1]
        )
    write_out(out, file_seqs)
    del seqs
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
    return ()


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

    Returns
    -------
    TYPE
        Description
    """
    fh = open(tmp_file, "w")
    fh.close()
    edge, ind = "", 0
    for id1 in edges:
        ida = id1
        if id1 in inverse:
            ida = inverse[id1]
            if ida.split("__")[0] in raw:
                ida = raw[ida.split("__")[0]]
        for id2 in edges[id1]:
            idb = id2
            if id2 in inverse:
                idb = inverse[id2]
                if idb.split("__")[0] in raw:
                    idb = raw[idb.split("__")[0]]
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
    del edges
    return ()


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

    Returns
    -------
    TYPE
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
    return ()


def reduce_edges(file_in, file_out):
    """Summary

    Parameters
    ----------
    file_in : TYPE
        Description
    file_out : TYPE
        Description

    Returns
    -------
    TYPE
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
    del done
    return ()


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

    Returns
    -------
    TYPE
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
    return ()


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

    Returns
    -------
    TYPE
        Description
    """
    (G, scale) = read_graphical_inputs(file_vertex, file_edges)
    output_cluster_file(G, cluster_file)
    return ()


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

    Returns
    -------
    TYPE
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
    return ()


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

    Returns
    -------
    TYPE
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
    return ()


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

    Returns
    -------
    TYPE
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
    return ()


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

    Returns
    -------
    TYPE
        Description
    """
    cluster_i(Sequence_file, tmp_file_1, edge_lengths)
    (s_sizes, cluster) = get_cluster_sizes_single(tmp_file_1)
    (seqs) = get_seqs_single(Sequence_file)
    get_similar_clusters(s_sizes, cluster, seqs, tmp_file + "_coclustered")
    (inv, coclust) = get_coclustered(tmp_file + "_coclustered")
    get_cluster_similarities_single(seqs, coclust, cluster, file_out, inv)
    return ()


##


def get_freq(id):
    """Summary

    Parameters
    ----------
    id : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    id = id.split("__")
    return int(id[1])


def intialise_files():
    """Initialise to output folder."""
    for d in [
        OUTFASTQ,
        OUTORTSEQ,
        OUTORTSEQTMP,
        OUTNET,
    ]:
        d.mkdir(exist_ok=True, parents=True)


def cram_to_bam(cram_path: Path, out_pre_qc_bam_path: Path):
    """Convert cram to bam.

    Parameters
    ----------
    cram_path : Path
        path to cram file.
    out_pre_qc_bam_path : Path
        path to output bam file.
    """
    cmd1 = [
        "samtools",
        "view",
        "-b",
        "-o",
        str(out_pre_qc_bam_path),
        str(cram_path),
    ]
    print(cmd)
    subprocess.run(cmd)


def bam_to_fastq(out_path: Path, source_path: Path, sample_id: str):
    """Summary

    Parameters
    ----------
    out_path : Path
        path to output folder.
    source_path : Path
        path of input file (cram or bam).
    sample_id : str
        name of sample. Also the prefix of the file.
    """
    if (
        len(
            glob(str(source_path.parent / "*.bam"))
            + glob(str(source.parent / "*.cram"))
        )
        != 0
    ):
        pre_qc_fastq1 = (
            out_path / "FASTQ_FILES" / f"Sequences_{sample_id}_1.fastq"
        )
        pre_qc_fastq2 = (
            out_path / "FASTQ_FILES" / f"Sequences_{sample_id}_2.fastq"
        )
        pre_qc_bam = out_path / "FASTQ_FILES" / f"Sequences_{sample_id}.bam"
        if len(glob(Path(source).parent / "*.cram")) != 0:
            cram_to_bam(cram_path=source, out_pre_qc_bam_path=pre_QC_bam)
        command1 = [
            "picard",
            "SamToFastq",
            f"I={str(pre_QC_bam)}",
            f"FASTQ={pre_QC_fastq1}",
            f"SECOND_END_FASTQ={pre_QC_fastq2}",
        ]
        subprocess.run(command1)


def qc_samples(
    out_path: Path,
    sample_id: str,
    min_length: Union[str, int] = 100,
    min_threshold: Union[str, int] = 32,
):
    """Perform QC on samples using QUASR and convert to fasta file with perl script.

    Parameters
    ----------
    out_path : Path
        path to output folder.
    sample_id : str
        name of sample. Also the prefix of the file.
    min_length : Union[str, int], optional
        minimum read length cutoff.
    min_threshold : Union[str, int], optional
        minimum median-read-quality cutoff.
    """
    reads1 = out_path / f"Sequences_{sample_id}_1.fastq"
    reads2 = out_path / f"Sequences_{sample_id}_2.fastq"
    quasr_qc_jar_path = EXTPATH / "QUASR_v7.01" / "qualityControl.jar"
    # see https://github.com/andrewjpage/QUASR for updated version
    cmd1 = [
        "java",
        "-jar",
        str(quasr_qc_jar_path),
        "-f",
        str(reads1),
        "-o",
        str(reads1.parent / reads1.stem),
        "-m",
        str(min_threshold),
        "-l",
        str(min_length),
    ]
    cmd2 = [
        "java",
        "-jar",
        str(quasr_qc_jar_path),
        "-f",
        str(reads2),
        "-o",
        str(reads2.parent / reads2.stem),
        "-m",
        str(min_threshold),
        "-l",
        str(min_length),
    ]
    cmd3 = [
        "perl",
        "-e",
        PERLCMD,
    ]
    subprocess.run(cmd1)
    subprocess.run(cmd2)
    subprocess.run(
        cmd3,
        stdin=open(reads1.with_suffix(".qc.fq"), "r"),
        stdout=open(reads1.with_suffix(".fasta"), "w"),
    )
    subprocess.run(
        cmd3,
        stdin=open(reads2.with_suffix(".qc.fq"), "r"),
        stdout=open(reads2.with_suffix(".fasta"), "w"),
    )


def prep_fastqs(out_path, source_path, sample_id, r1pattern, r2pattern):
    """Prepare fastqs for input into the script.

    Parameters
    ----------
    out_path : Path
        location of output folder.
    source_path : Path
        location of R1 of fastq file.
    sample_id : str
        name of sample. Also the prefix of the file.
    r1pattern : str
        suffix pattern before .fastq to try and match for R1
    r2pattern : str
        suffix pattern before .fastq to try and match for R2
    """
    if re.search(r1pattern, str(source_path)):
        r1_original = Path(source_path)
        r2_original = Path(re.sub(r1pattern, r2pattern, str(source_path)))
    else:
        raise ValueError(
            "Your input file {} does not contain the {} pattern.".format(
                str(source_path), r1pattern
            )
        )
    copy_prepped_fastq(
        out_path=out_path,
        fastq_path=r1_original,
        sample_id=sample_id,
        read_num="1",
    )
    copy_prepped_fastq(
        out_path=out_path,
        fastq_path=r2_original,
        sample_id=sample_id,
        read_num="2",
    )


def copy_prepped_fastq(
    out_path: Path, fastq_path: Path, sample_id: str, read_num: str
):
    """Copy and unzip fastq files.

    Parameters
    ----------
    out_path : Path
        location of output folder
    fastq_path : Path
        location of original fastq file.
    sample_id : str
        name of sample. Also the prefix of the file.
    read_num : str
        `1` or `2` to specify read 1 or read 2.
    """
    if fastq_path.suffix == ".gz":
        extension = ".gz"
    else:
        extension = ""
    new_fastq = (
        out_path
        / "FASTQ_FILES"
        / f"Sequences_{sample_id}_{read_num}.fastq{extension}"
    )
    shutil.copy(fastq_path, new_fastq)
    if fastq_path.suffix == ".gz":
        cmdg = ["gunzip", "-f", f"{new_fastq}"]
        subprocess.run(cmdg)


###########################
out_path = Path(sys.argv[1])
sample_id = sys.argv[2]
barcode_group = sys.argv[3]
gene = sys.argv[4]
paired = sys.argv[5]
species = sys.argv[6]
source = sys.argv[7]
length = sys.argv[8]
primer_file = sys.argv[9]
method = sys.argv[10]
command_source = sys.argv[11]
command_source = command_source.split(",")
if len(sys.argv) > 13:
    reverse_primer_group = sys.argv[13]
else:
    reverse_primer_group = "OTHER"
print("Reverse primer group: ", reverse_primer_group)

OUTFASTQ = out_path / "FASTQ_FILES"
OUTORTSEQ = out_path / "ORIENTATED_SEQUENCES"
OUTORTSEQTMP = OUTORTSEQ / "TMP"
OUTNET = OUTORTSEQ / "NETWORKS"

# Files for QC and filtering
Seq_file1 = OUTFASTQ / f"Sequences_{sample_id}_1.fasta"
Seq_file2 = OUTFASTQ / f"Sequences_{sample_id}_2.fasta"
tmp_Tmp_file = OUTORTSEQTMP / f"Untrimmed_{sample_id}.fasta"
trim1 = OUTORTSEQTMP / f"trimmed_orientated_all_{sample_id}.fasta"
trim2 = OUTORTSEQTMP / f"Filtered_J_{sample_id}.fasta"
trim3 = OUTORTSEQTMP / f"Filtered_reduced_{sample_id}.fasta"
Fail_file = OUTFASTQ / f"Fail_filtered_{sample_id}.fasta"
primer_tag_file = (
    OUTORTSEQTMP / f"Barcode_filtering_information_{sample_id}.txt"
)
primer_tag_file_count = OUTORTSEQTMP / f"All_barcodes_{sample_id}.txt"
Filtered_out1 = OUTORTSEQ / f"Filtered_ORFs_sequences_all_{sample_id}.fasta"
nn_orf_filtered = (
    OUTORTSEQ / f"Nucleotsample_ide_ORF_filtered_all_{sample_id}.fasta"
)
tmp_file_orf = OUTORTSEQTMP / f"Blast_matching_{sample_id}"
filtering_report = OUTORTSEQ / f"Filtering_report_{sample_id}.txt"
# Files for clustering
att_file = OUTNET / f"Vertex_relations_{sample_id}.txt"
file_vertex = OUTNET / f"Att_{sample_id}.txt"
file_edges = OUTNET / f"Edges_{sample_id}.txt"
cluster_file = OUTNET / f"Cluster_sample_identities_{sample_id}.txt"
Reduced_file = OUTNET / f"Fully_reduced_{sample_id}.fasta"
checked_edges = OUTNET / f"Checked_edges_{sample_id}.txt"
plot_sample_ids_file = OUTNET / f"Plot_sample_ids_{sample_id}.txt"
file_seqs = OUTNET / f"Sequences_{sample_id}.txt"
tmp_file0 = OUTNET / f"Decon_0_{sample_id}.txt"
tmp_pre = OUTNET / "Pre_tmp_{sample_id}"
tmp_file_1 = OUTNET / f"NN_Tmp_cluster_{sample_id}.1"

tmp_file = OUTNET / f"NN_Tmp_cluster_{sample_id}."

# Reference files
refv = str(LIBPATH / f"Reference_nn_{species}_{gene}V.fasta")
refj = str(LIBPATH / f"Reference_nn_{species}_{gene}J.fasta")
refvp = str(LIBPATH / f"Reference_protein_{species}_{gene}V.fasta")
refjp = str(LIBPATH / f"Reference_protein_{species}_{gene}J.fasta")
ref_const = str(LIBPATH / f"Reference_nn_{species}_{gene}_constant_exon1.fasta")

# change here if necessary
R1PATTERN = "_R1_001"
R2PATTERN = "_R2_001"

# Commands
if command_source.count("1") != 0:
    intialise_files()
    if (
        len(
            glob(str(Path(source).parent / "*.bam*"))
            + glob(str(Path(source).parent / "*.cram"))
        )
        != 0
    ):
        bam_to_fastq(out_path=out_path, source_path=source, sample_id=sample_id)
    elif (
        len(
            glob(str(Path(source).parent / "*.fastq"))
            + glob(str(Path(source).parent / "*.fastq.gz"))
            + glob(str(Path(source).parent / "*.fq*"))
            + glob(str(Path(source).parent / "*.fq.gz"))
        )
        != 0
    ):
        # rename them to Seqeuence_{sample_id}_1.fastq Seqeuence_{sample_id}_2.fastq
        prep_fastqs(
            out_path=out_path,
            source_path=source,
            sample_id=sample_id,
            r1pattern=R1PATTERN,
            r2pattern=R2PATTERN,
        )
    qc_samples(
        out_path=OUTFASTQ,
        sample_id=sample_id,
        min_length=length,
        min_threshold=MIN_QUAL,
    )


# Tip: it is good to check all the fasta files in the FASTQ_FILES directory have
# been made correctly at this point (with non-zero number of lines)

# Filtering and processing reads
if command_source.count("2") != 0:
    if gene.count("IG") != 0:
        get_paired_reads_overlapping(
            Seq_file1, Seq_file2, tmp_Tmp_file, gene, paired, sample_id, method
        )
        trim_sequences_bcr_tcr(
            tmp_Tmp_file,
            Fail_file,
            trim1,
            gene,
            paired,
            species,
            primer_file,
            primer_tag_file,
            tmp_file,
            primer_tag_file_count,
            sample_id,
            ref_const,
            reverse_primer_group,
        )
        filter_igj_genes(
            trim1,
            trim2,
            refj,
            primer_file,
            ref_const,
            primer_tag_file_count,
        )
        reduce_sequences(trim2, trim3, primer_file)
        orf_calculation_single(
            trim3,
            Filtered_out1,
            nn_orf_filtered,
            out_path,
            gene,
            refv,
            refj,
            refvp,
            refjp,
            tmp_file_orf,
        )
        get_read_report(
            Seq_file1,
            Seq_file2,
            tmp_Tmp_file,
            trim1,
            nn_orf_filtered,
            filtering_report,
            sample_id,
            species,
            gene,
            out_path,
            primer_tag_file_count,
            primer_file,
            method,
            barcode_group,
        )

# Clustering reads
if command_source.count("3") != 0:
    generate_networks(
        nn_orf_filtered, tmp_file_1, EDGE_LENGTHS, tmp_file, att_file
    )
    deconvolute_edges(
        nn_orf_filtered,
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
