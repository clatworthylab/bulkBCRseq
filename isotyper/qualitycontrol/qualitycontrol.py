#!/usr/bin/env python


def get_paired_reads_overlapping(file1, file2, outfile, paired):
    """Summary

    Parameters
    ----------
    file1 : TYPE
        Description
    file2 : TYPE
        Description
    outfile : TYPE
        Description
    paired : TYPE
        Description
    """
    seqs1 = get_sequences(file1)
    seqs2 = get_sequences(file2)
    print("Forward reads:", len(seqs1), "Reverse reads:", len(seqs2))
    out, ind = "", 0
    fh = open(outfile, "w")
    fh.close()
    tot, f_joining, total, length, fail_no_pair = 0, 0, 0, 30, 0
    for seq_id in seqs1:
        if seq_id in seqs2:
            total = total + 1
            (seq, failed) = join_reads(seqs1[seq_id], seqs2[seq_id], length)
            if failed == 0 and len(seq) > 180:
                ind, tot = ind + 1, tot + 1
                out = out + ">" + seq_id + "\n" + seq + "\n"
            else:
                seq2 = reverse_comp(seqs2[seq_id])
                (seq, failed) = join_reads(seqs1[seq_id], seq2, length)
                if failed == 0 and len(seq) > 180:
                    ind, tot = ind + 1, tot + 1
                    out = out + ">" + seq_id + "\n" + seq + "\n"
                else:
                    f_joining = f_joining + 1
            if ind > 100:
                write_out(out, outfile)
                out, ind = "", 0
        else:
            fail_no_pair = fail_no_pair + 1
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


def read_untrimmed_file_single(
    tmp_Tmp_file,
    Fail_file,
    Output_trim,
    paired,
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
    paired : TYPE
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
    """
    for f in [primer_tag_file_count]:
        fh = open(f, "w")
        fh.close()
    fh = open(tmp_Tmp_file, "r")
    minl, maxl = 120, 1000000
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
    paired,
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
    """
    read_untrimmed_file_single(
        tmp_Tmp_file,
        Fail_file,
        Output_trim,
        paired,
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
        h = (
            h.split(":")[0].split(READ_NUMBER_DIVISION)[0]
            + READ_NUMBER_DIVISION
            + str(sum(u_freq))
        )
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
            ind = ind + 1
            if ind > 100:
                write_out(out, Output_trim)
                out, ind = "", 0
    write_out(out, Output_trim)


def trim_sequences_bcr_tcr(
    tmp_Tmp_file,
    Fail_file,
    Output_trim,
    paired,
    primer_tag_file,
    tmp_file,
    primer_tag_file_count,
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
    paired : TYPE
        Description
    primer_tag_file : TYPE
        Description
    tmp_file : TYPE
        Description
    primer_tag_file_count : TYPE
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
        PRIMER_FILE
    )
    fh_out = open(Output_trim, "w")
    fh_out.close()
    fh_out = open(Fail_file, "w")
    fh_out.close()
    if barcoded_j == 1 and barcoded_v == 0:
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
            paired,
            primer_tag_file_count,
            ref_const,
            v_ref,
        )


get_paired_reads_overlapping(Seq_file1, Seq_file2, tmp_Tmp_file, paired)
trim_sequences_bcr_tcr(
    tmp_Tmp_file,
    Fail_file,
    trim1,
    paired,
    primer_tag_file,
    tmp_file,
    primer_tag_file_count,
    ref_const,
    reverse_primer_group,
)
filter_igj_genes(
    trim1,
    trim2,
    refj,
    ref_const,
    primer_tag_file_count,
)
reduce_sequences(trim2, trim3)
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
    barcode_group,
)
