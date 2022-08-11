#!/usr/bin/env python
import subprocess

from operator import add
from pathlib import Path
from typing import Dict, Tuple, List

from isotyper.utilities._utils import (
    check_fasta_not_empty,
    cluster_i,
    fasta_iterator,
    get_codons,
    get_freq,
    get_match,
    join_reads,
    reverse_comp,
    translate,
    Tree,
    write_out,
)
from isotyper.utilities._settings import (
    READ_NUMBER_DIVISION,
    THRESHOLD_BARCODE,
    THRESHOLD_GENE_SCORE,
)
from isotyper.qualitycontrol._settings import *


def get_sequences(file: Path) -> Dict:
    """Extract sequences from fasta file.

    Parameters
    ----------
    file : Path
        path to fasta file.

    Returns
    -------
    Dict
        dictionary holding fasta header and sequence as key and record respectively.
    """
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        header = header.replace(":", "").split()[0].replace("-", "")
        header = header.split("#")[0].split("/")[0]
        seqs[header] = sequence
    fh.close()
    return seqs


def get_paired_reads_overlapping(
    seq_file1: Path, seq_file2: Path, outfile: Path
):
    """Summary

    Parameters
    ----------
    seq_file1 : Path
        path to pre-QC read 1 in fasta format.
    seq_file2 : Path
        path to pre-QC read 2 in fasta format.
    outfile : Path
        path to untrimmed seqeunces in fasta format.
    """
    seqs1 = get_sequences(seq_file1)
    seqs2 = get_sequences(seq_file2)
    print("Forward reads:", len(seqs1), "Reverse reads:", len(seqs2))
    out, ind = "", 0
    fh = open(outfile, "w")
    fh.close()
    tot, f_joining, total, length, fail_no_pair = 0, 0, 0, 30, 0
    for seq_id in seqs1:
        if seq_id in seqs2:
            total += 1
            (seq, failed) = join_reads(seqs1[seq_id], seqs2[seq_id], length)
            if failed == 0 and len(seq) > 180:
                ind += 1
                tot += 1
                out = out + ">" + seq_id + "\n" + seq + "\n"
            else:
                seq2 = reverse_comp(seqs2[seq_id])
                (seq, failed) = join_reads(seqs1[seq_id], seq2, length)
                if failed == 0 and len(seq) > 180:
                    ind += 1
                    tot += 1
                    out = out + ">" + seq_id + "\n" + seq + "\n"
                else:
                    f_joining += 1
            if ind > 100:
                write_out(out, outfile)
                out, ind = "", 0
        else:
            fail_no_pair += 1
    write_out(out, outfile)
    del out, ind
    print(
        SAMPLE_ID + "\tTotal sequences:",
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


def trim_sequences_bcr_tcr(
    untrimmed_file: Path,
    output_trim: Path,
    primer_tag_file: Path,
    primer_tag_file_count: Path,
    ref_const: Path,
):
    """Summary

    Parameters
    ----------
    untrimmed_file : Path
        path to untrimmed fasta file.
    output_trim : Path
        path to output fasta file.
    primer_tag_file : Path
        path to primer tag table.
    primer_tag_file_count : Path
        path to another primer tag table.
    ref_const : Path
        path to constant reference file.
    """
    forward, reverse, barcoded_j, barcoded_v, v_ref = get_primers_split(
        PRIMER_FILE
    )
    fh_out = open(output_trim, "w")
    fh_out.close()
    fh_out = open(FAIL_FILE, "w")
    fh_out.close()
    if barcoded_j == 1 and barcoded_v == 0:
        single_j_barcoded_trimming_clustered(
            forward=forward,
            reverse=reverse,
            barcoded_j=barcoded_j,
            barcoded_v=barcoded_v,
            untrimmed_file=untrimmed_file,
            output_trim=output_trim,
            primer_tag_file=primer_tag_file,
            primer_tag_file_count=primer_tag_file_count,
            ref_const=ref_const,
            v_ref=v_ref,
        )


def single_j_barcoded_trimming_clustered(
    forward,
    reverse,
    barcoded_j,
    barcoded_v,
    untrimmed_file: Path,
    output_trim: Path,
    primer_tag_file: Path,
    primer_tag_file_count: Path,
    ref_const: Path,
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
    untrimmed_file : Path
        Description
    output_trim : Path
        Description
    primer_tag_file : Path
        Description
    primer_tag_file_count : Path
        Description
    ref_const : TYPE
        Description
    v_ref : TYPE
        Description
    """
    read_untrimmed_file_single(
        untrimmed_file=untrimmed_file,
        primer_tag_file_count=primer_tag_file_count,
        forward=forward,
        reverse=reverse,
        v_ref=v_ref,
    )
    check_barcodes_malbac(
        primer_tag_file_count=primer_tag_file_count,
        primer_tag_file=primer_tag_file,
        output_trim=output_trim,
    )
    separate_sequences(
        primer_tag_file_count=primer_tag_file_count,
        primer_tag_file=primer_tag_file,
        output_trim=output_trim,
        ref_const=ref_const,
    )


def check_barcodes_malbac(
    primer_tag_file_count: Path, primer_tag_file: Path, output_trim: Path
):
    """Summary

    Parameters
    ----------
    primer_tag_file_count : TYPE
        Description
    primer_tag_file : TYPE
        Description
    output_trim : TYPE
        Description
    """
    header_txt = (
        "#ID\tnumber_of_reads\ttotal_reads_with_BC\tJ_barcode\tV_barcode\t"
        + "bp_mismatches_to_consensus\tBC_accepted\tconsensus\tsequence\n"
    )
    outs = ["", header_txt]
    files = [output_trim, primer_tag_file]
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
                f += int(h.split(":")[1])
            u_seq, u_freq, u_header = (
                u_seq + [s],
                u_freq + [f],
                u_header + [h],
            )
            total_tags += f
            ind += 1
        f = sum(u_freq)
        h = (
            h.split(":")[0].split(READ_NUMBER_DIVISION)[0]
            + READ_NUMBER_DIVISION
            + str(sum(u_freq))
        )
        if sum(u_freq) > 20:
            print("BC group size:\t", len(u_freq), t1.split())
        if len(u_freq) == 1:
            passed_seqs_total += f
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
                passed_seqs_total += f
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
                    len(u_seq[i]),
                    ids,
                    sum(u_freq),
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
                    passed_seqs_total += f
                else:
                    fail_less_than_threshold += 1
            else:
                consensus, pass_consensus = get_consensus_sequence_cluster(
                    u_seq, u_freq
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
                    passed_seqs_total += f
                else:
                    fail_less_than_threshold += 1
        if ind > 200:
            write_out(outp, primer_tag_file)
            outp, ind = "", 0
    write_out(outp, primer_tag_file)
    outp, ind = "", 0
    print(total_tags, passed_seqs_total, fail_less_than_threshold)


def get_consensus_sequence_large(out_cluster, l1, ids, sum_u_freq):
    """Summary

    Parameters
    ----------
    out_cluster : TYPE
        Description
    l1 : TYPE
        Description
    ids : TYPE
        Description
    sum_u_freq : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    fh = open(FAIL_FILE, "w")
    fh.write(out_cluster)
    fh.close()
    cluster_i(
        input_file=FAIL_FILE,
        clustered_file=FAIL_FILE.with_suffix(FAIL_FILE.suffix + "cls"),
        identity=(l1 - 3.0) / l1,
    )
    cluster = Tree()
    fh = open(FAIL_FILE.with_suffix(FAIL_FILE.suffix + "cls.bak.clstr"), "r")
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
        consensus = get_consensus_sequence(seq_align, s_freq)
        if consensus.count("_") == 0 and len(consensus) > 3:
            pass_consensus = 1
        else:
            pass_consensus = 0
    return (consensus, pass_consensus)


def get_consensus_sequence(u_seq, u_freq):
    """Summary

    Parameters
    ----------
    u_seq : TYPE
        Description
    u_freq : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    out = ""
    for i in range(0, len(u_seq)):
        out = out + ">" + str(i) + "\n" + u_seq[i] + "\n"
    fh = open(TMP_FILE, "w")
    fh.write(out)
    fh.close()
    main_cmd = [
        "mafft",
        "--quiet",
        "--auto",
        "--retree",
        "2",
    ]
    if len(u_seq) > 2000:
        main_cmd = main_cmd + ["--parttree"]
    main_cmd = main_cmd + [f"{TMP_FILE}"]
    subprocess.run(
        main_cmd,
        stdout=open(TMP_FILE.with_suffix(".aligned"), "w"),
    )
    max_seqs = {}
    pfh = check_fasta_not_empty(TMP_FILE.with_suffix(".aligned"))
    if pfh == 1:
        fh = open(TMP_FILE.with_suffix(".aligned"), "r")
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
    else:
        consensus = ""
    return consensus


# TODO: this function looks a lot that is duplicated in get_consensus_sequence. maybe can combine?
def get_consensus_sequence_cluster(u_seq, u_freq):
    """Summary

    Parameters
    ----------
    u_seq : TYPE
        Description
    u_freq : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    out = ""
    for i in range(0, len(u_seq)):
        out = out + ">" + str(i) + "\n" + u_seq[i] + "\n"
    fh = open(TMP_FILE, "w")
    fh.write(out)
    fh.close()
    main_cmd = [
        "mafft",
        "--quiet",
        "--auto",
        "--retree",
        "2",
    ]
    if len(u_seq) > 2000:
        main_cmd = main_cmd + ["--parttree"]
    main_cmd = main_cmd + [f"{TMP_FILE}"]
    subprocess.run(
        main_cmd,
        stdout=open(TMP_FILE.with_suffix(".aligned"), "w"),
    )
    fh = open(TMP_FILE.with_suffix(".aligned"), "r")
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


def read_untrimmed_file_single(
    untrimmed_file: Path,
    primer_tag_file_count: Path,
    forward,
    reverse,
    v_ref,
):
    """Summary

    Parameters
    ----------
    untrimmed_file : Path
        Description
    primer_tag_file_count : Path
        Description
    forward : TYPE
        Description
    reverse : TYPE
        Description
    v_ref : TYPE
        Description
    """
    fh = open(primer_tag_file_count, "w")
    fh.close()
    fh = open(untrimmed_file, "r")
    minl, maxl = 120, 1000000
    seqs1, t = Tree(), 0
    for header, seq in fasta_iterator(fh):
        seq = seq.upper()
        seqs1[seq][header].value = 1
        t += 1
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
        total += 1
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
            pass_r += 1
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
                                p += 1
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
            pass_f += 1
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
                    pass_all += 1
                    ind += 1
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


def separate_sequences(
    primer_tag_file_count: Path,
    primer_tag_file: Path,
    output_trim: Path,
    ref_const,
):
    """Summary

    Parameters
    ----------
    primer_tag_file_count : Path
        Description
    primer_tag_file : Path
        Description
    output_trim : Path
        Description
    ref_const : TYPE
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
                id = l[0].split(READ_NUMBER_DIVISION)[0]
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
                scores[0][0] > THRESHOLD_GENE_SCORE
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
                + id.split(READ_NUMBER_DIVISION)[0]
                + READ_NUMBER_DIVISION
                + "_".join(map(str, f))
                + "|"
                + header
                + "\n"
                + s
                + "\n"
            )
            ind += 1
            if ind > 100:
                write_out(out, output_trim)
                out, ind = "", 0
    write_out(out, output_trim)


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


def orf_calculation_single(
    output_trim: Path,
    filtered_out: Path,
    nn_orf_filtered: Path,
    ref_protein: Path,
):
    """get protein and nucleotide sequences

    Parameters
    ----------
    output_trim : Path
        Description
    filtered_out : Path
        Description
    nn_orf_filtered : Path
        Description
    ref_protein : Path
        Description
    """
    get_protein_sequences(
        output_trim=output_trim,
        filtered_out=filtered_out,
        nn_orf_filtered=nn_orf_filtered,
        ref_protein=ref_protein,
    )
    get_nucleotide_sequences(
        output_trim=output_trim,
        filtered_out=filtered_out,
        nn_orf_filtered=nn_orf_filtered,
    )


def get_sequences_ref(file: Path, word: str) -> Dict:
    """Summary

    Parameters
    ----------
    file : Path
        Description
    word : str
        Description

    Returns
    -------
    Dict
        Description
    """
    fh = open(file, "r")
    dt = Tree()
    for header, sequence in fasta_iterator(fh):
        for i in range(1, len(sequence) - word - 1):
            dt[sequence[i : (i + word)]][header].value = 1
    fh.close()
    return dt


def calculate_orf_length(codon: Dict, sequence: str) -> Tuple[List, int]:
    """Summary

    Parameters
    ----------
    codon : Dict
        dictionary of codons.
    sequence : str
        sequence to translate.

    Returns
    -------
    Tuple[List, int]
        list of accepted orfs.
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


def get_max_orf(orf_all: str, word: str, dt: Dict) -> Tuple:
    """Summary

    Parameters
    ----------
    orf_all : str
        Description
    word : str
        Description
    dt : Dict
        Description

    Returns
    -------
    Tuple
        Description
    """
    score = []
    codon = []
    if len(orf_all) > 1:
        for si in orf_all:
            s = si[0]
            sc = 0
            for i in range(1, len(s) - word):
                if s[i : (i + word)] in dt:
                    sc += 1
            score.append(sc)
            codon.append(si[1])
        maxi = score.index(max(score))
        orf = orf_all[maxi][0]
        mscore = max(score)
        codon = orf_all[maxi][1]
    else:
        for si in orf_all:
            orf = si[0]
            mscore = 1000
            codon = si[1]
    return (orf, mscore, codon)


def get_protein_sequences(
    output_trim: Path,
    filtered_out: Path,
    nn_orf_filtered: Path,
    ref_protein: Path,
):
    """write protein sequences.

    Parameters
    ----------
    output_trim : Path
        output reduced fasta file.
    filtered_out : Path
        filtered fasta file
    nn_orf_filtered : Path
        filtered nucleotide fasta file.
    ref_protein : Path
        reference fasta file.
    """
    fh = open(filtered_out, "w")
    fh.close()
    fh = open(nn_orf_filtered, "w")
    fh.close()
    codon = get_codons(codon_file=CODON_FILE)
    (word) = 4
    (dict1) = get_sequences_ref(file=ref_protein, word=word)
    fh = open(output_trim, "r")
    number_seqs, number_accepted, batch, indb = 0, 0, 1000, 0
    # V_region, seqs, ORFa, ORFb = '', Tree(), {}, {}
    seqs, ORFa = Tree(), {}
    batch_number, orf_fail, p_seq = 0, 0, {}
    out = ""
    for header, seq in fasta_iterator(fh):
        if seq.count("N") == 0:
            seq = seq.replace("-", "")
            header = header.split("#")[0]
            if header.count(READ_NUMBER_DIVISION) == 0:
                freq = get_freq(header)
                header = (
                    header.replace(":", "") + READ_NUMBER_DIVISION + str(freq)
                )
            (ORF1, accept1) = calculate_orf_length(codon=codon, sequence=seq)
            ORFa[header] = ORF1
            number_seqs += 1  # freq
            if accept1 == 0:
                orf_fail += 1
            if accept1 != 0:
                if len(ORF1) > 1:
                    min_score, found = 2, 0
                    (O1, score1, codon1) = get_max_orf(
                        orf_all=ORF1, word=word, dt=dict1
                    )
                    if score1 > min_score:
                        min_score = score1
                        p_seq[header] = O1
                        seq = seq[codon1 : len(seq)]
                        found = 1
                    if found == 1:
                        number_accepted += 1
                else:
                    p_seq[header] = ORF1[0][0]
                    number_accepted += 1
                    seq = seq[ORF1[0][1] : len(seq)]
                seqs[seq][header].value = 1
                indb += 1
                if indb >= batch:
                    batch_number += 1
                    out = ""
                    for seq_id in p_seq:
                        out = out + ">" + seq_id + "\n" + p_seq[seq_id] + "\n"
                    write_out(out, filtered_out)
                    p_seq = {}
                    out = ""
    fh.close()
    write_out(out, filtered_out)
    out, ind = "", 0
    for seq_id in p_seq:
        out = out + ">" + seq_id + "\n" + p_seq[seq_id] + "\n"
        ind += 1
        if ind > 500:
            write_out(out, filtered_out)
            out, ind = "", 0
    write_out(out, filtered_out)


def get_nucleotide_sequences(
    output_trim,
    filtered_out,
    nn_orf_filtered,
):
    """write nucleotide sequences.

    Parameters
    ----------
    output_trim : Path
        output reduced fasta file.
    filtered_out : Path
        filtered fasta file
    nn_orf_filtered : Path
        filtered nucleotide fasta file.
    """
    fh = open(nn_orf_filtered, "w")
    fh.close()
    fh = open(filtered_out, "r")
    ids = {}
    for header, sequence in fasta_iterator(fh):
        ids[header] = 1
    fh.close()
    out, ind = "", 0
    fh = open(output_trim, "r")
    done = {}
    for header, sequence in fasta_iterator(fh):
        if header in ids:
            out = out + ">" + header + "\n" + sequence + "\n"
            done[header] = 1
            ind += 1
            if ind > 500:
                write_out(out, nn_orf_filtered)
                out, ind = "", 0
    write_out(out, nn_orf_filtered)
    fh.close()
    print(len(ids), len(done))


def filter_igj_genes(trim1: Path, trim2: Path, refj):
    """Filter IGHJ genes based on quality

    Parameters
    ----------
    trim1 : Path
        path to all trimmed sequences.
    trim2 : Path
        path to all trimmed sequences after blast for J gene.
    refj : TYPE
        IGHJ reference.
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
        inf = 0
        out = out + ">" + header + "\n" + sequence[inf : len(sequence)] + "\n"
        seqs[header] = sequence
        ind += 1
        batch += 1
        c += 1
        if batch >= batch_size:
            blast_match_j(
                out=out,
                seqs=seqs,
                trim1=trim1,
                trim2=trim2,
                refj=refj,
                e_value=e_value,
            )
            out, batch, seqs = "", 0, {}
    fh.close()
    if len(out) > 2:
        blast_match_j(
            out=out,
            seqs=seqs,
            trim1=trim1,
            trim2=trim2,
            refj=refj,
            e_value=e_value,
        )
        out, batch = "", 0


def blast_match_j(
    out: str, seqs: Dict, trim1: Path, trim2: Path, refj: Path, e_value: float
):
    """Use blast to QC J gene assignments.

    Parameters
    ----------
    out : str
        Description
    seqs : Dict
        Description
    trim1 : Path
        path to all trimmed sequences.
    trim2 : Path
        path to all trimmed sequences after blast for J gene.
    refj : Path
        path to j reference file.
    e_value : float
        e-value cut off.
    """
    blasted_j_file = trim1.with_suffix(trim1.suffix + "_blast_J")
    blasted_j_results = trim1.with_suffix(trim1.suffix + "_blast_J_results")
    fh = open(blasted_j_file, "w")
    fh.write(out)
    fh.close()
    command1 = [
        "blastn",
        "-num_threads",
        "10",
        "-db",
        f"{str(refj)}",
        "-evalue",
        f"{str(e_value)}",
        "-query",
        f"{str(blasted_j_file)}",
        "-out",
        f"{str(blasted_j_results)}",
        "-num_alignments",
        "1",
        "-outfmt",
        "6",
        "-word_size",
        "4",
    ]
    subprocess.run(command1)
    fh = open(blasted_j_results, "r")
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


def reduce_sequences(trim2, trim3):
    """Summary

    Parameters
    ----------
    trim2 : TYPE
        Description
    trim3 : TYPE
        Description
    """
    minl = 185  # change for shorter runs
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
        for seq_id in seqs[seq]:
            f1 = list(
                map(
                    int,
                    seq_id.split(READ_NUMBER_DIVISION)[1]
                    .split("|")[0]
                    .split("_"),
                )
            )
            f = list(map(add, f, f1))
        header = (
            ">"
            + seq_id.split(READ_NUMBER_DIVISION)[0]
            + READ_NUMBER_DIVISION
            + "_".join(map(str, f))
            + "|"
            + head
        )
        out = out + header + "\n" + seq + "\n"
        ind += 1
        if ind > 500:
            write_out(out, trim3)
            out, ind = "", 0
    write_out(out, trim3)


def count_barcodes(primer_tag_file: Path) -> Tuple[int, int, int]:
    """Count number of barcodes and unique sequences.

    Parameters
    ----------
    primer_tag_file : Path
        path to barcode filtering information table.

    Returns
    -------
    Tuple[int, int, int]
        number of barcodes and unique sequences.
    """
    fh = open(primer_tag_file, "r")
    bc_uniq, bc_count, bc_failed = {}, 0, 0
    uniq_seq = 0
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            bc = l[3] + ":" + l[4]
            if l[1] != "NA":
                freq = int(l[1])  # int(l[0].split(READ_NUMBER_DIVISION)[1])
                uniq_seq += 1
                bc_uniq[bc] = 1
                bc_count += freq
                if l[6] != "YES":
                    bc_failed += 1
    fh.close()
    n_barcodes, uniq_sequences, total_sequences_included_before_bc = (
        len(bc_uniq),
        uniq_seq,
        bc_count,
    )
    return (n_barcodes, uniq_sequences, total_sequences_included_before_bc)


def get_number_sequences(file: Path) -> int:
    """Count number of sequences (umi) in file

    Parameters
    ----------
    file : Path
        Input fasta file.

    Returns
    -------
    int
        number of sequences.
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
            count += 1
        fh.close()
    else:
        count = -1
    return count


def get_reduced_number_sequences_multi_constants(file: Path) -> int:
    """Count number of sequences (umi) with multiple contant gene assignments.

    Parameters
    ----------
    file : Path
        Input fasta file.

    Returns
    -------
    int
        number of sequences.
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
                map(
                    int,
                    header.split(READ_NUMBER_DIVISION)[1]
                    .split("|")[0]
                    .split("_"),
                )
            )
        fh.close()
    return count


def get_read_report(
    seq_file1: Path,
    seq_file2: Path,
    untrimmed_file: Path,
    trim1: Path,
    nn_orf_filtered: Path,
    filtering_report: Path,
    sample_id: str,
    species: str,
    out_path: Path,
):
    """Generate report.

    Parameters
    ----------
    seq_file1 : Path
        pre qc read 1 fasta.
    seq_file2 : Path
        pre qc read 2 fasta.
    untrimmed_file : Path
        untrimmed paired fasta
    trim1 : Path
        trimmed fasta
    nn_orf_filtered : Path
        filtered fasta
    filtering_report : Path
        filtering report
    sample_id : str
        sample name
    species : str
        organism
    out_path : Path
        location of output.
    """
    barcoded_j, barcoded_v = 1, 0
    print(barcoded_j, barcoded_v)
    (
        n_barcodes,
        uniq_sequences,
        total_sequences_included_before_bc,
    ) = count_barcodes(primer_tag_file=primer_tag_file)
    print(n_barcodes, uniq_sequences, total_sequences_included_before_bc)
    (
        raw1,
        raw2,
        joined,
        number_unique_seqs,
        gene_matching,
        orf,
        number_unique_seqs,
    ) = (
        get_number_sequences(file=seq_file1),
        get_number_sequences(file=seq_file2),
        get_number_sequences(file=untrimmed_file),
        get_number_sequences(file=nn_orf_filtered),
        get_reduced_number_sequences_multi_constants(file=trim1),
        get_reduced_number_sequences_multi_constants(file=nn_orf_filtered),
    )
    out = (
        "Directory\tSample\tSpecies\tGene\t% reads retained\tN raw reads (1)\t"
        + "N raw reads (2)\tN joined reads\tN reads with BCs\tN uniq BCs\tN reads gene matched\t"
        + "N reads w/t ORF\tUnique sequences\n"
    )
    orf_perc = orf * 100.0 / min([raw1, raw2])
    if min([raw1, raw2]) == -1:
        orf_perc = "NA"
    out = (
        out
        + str(out_path)
        + "\t"
        + sample_id
        + "\t"
        + species
        + "\t"
        + "IGH"  # always IGH in Clatworthy Lab
        + "\t"
        + str(orf_perc)
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


def main():
    """main function in step 2."""
    get_paired_reads_overlapping(
        seq_file1=SEQ_FASTA_FILE1,
        seq_file2=SEQ_FASTA_FILE2,
        outfile=UNTRIMMED_FASTA,
    )
    trim_sequences_bcr_tcr(
        untrimmed_file=UNTRIMMED_FASTA,
        output_trim=TRIM1_ALL,
        primer_tag_file=PRIMER_TAG_FILE,
        primer_tag_file_count=PRIMER_TAG_FILE_COUNT,
        ref_const=REF_CONST,
    )
    filter_igj_genes(
        trim1=TRIM1_ALL,
        trim2=TRIM2_J,
        refj=REFJ,
    )
    reduce_sequences(trim2=TRIM2_J, trim3=TRIM3_RED)
    orf_calculation_single(
        output_trim=TRIM3_RED,
        filtered_out=FILTERED_OUT,
        nn_orf_filtered=FILTERED_OUT_NT,
        ref_protein=REFVP,
    )
    get_read_report(
        seq_file1=SEQ_FASTA_FILE1,
        seq_file2=SEQ_FASTA_FILE2,
        untrimmed_file=UNTRIMMED_FASTA,
        trim1=TRIM1_ALL,
        nn_orf_filtered=FILTERED_OUT_NT,
        filtering_report=FILTERING_REPORT,
        sample_id=SAMPLE_ID,
        species=ORG,
        out_path=OUT_PATH,
    )


if __name__ == "__main__":
    main()
