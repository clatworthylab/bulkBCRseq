def single_j_barcoded_trimming(regions_J, J_primer, J1, J2, V_primer, V1, V2,
                               rc, Tmp_file, Fail_file, Output_trim, gene,
                               paired, species, primer_file, primer_tag_file,
                               tmp_file, sample, primer_tag_file_count,
                               threshold_barcode, vidprimer, jidprimer, inside,
                               ref_const, universal_rev, reverse_primer_group):
    (js, vs) = (len(J_primer), len(V_primer))
    minl, maxl = 150, 1000000
    if (gene == "HEAVY" or gene == "IGH"):
        minl = 150
    if (gene == "KAPPA" or gene == "IGK"):
        minl = 110
    if (gene == "LAMBDA" or gene == "IGL"):
        (minl, maxl) = (90, 150)
    # J_found, v_found, tot, indexing = 0, 0, 0, 0
    print(primer_file, inside)
    read_untrimmed_file(Tmp_file, J_primer, rc, J1, J2, regions_J, sample,
                        maxl, minl, js, vs, V_primer, V1, V2,
                        primer_tag_file_count, Fail_file, vidprimer, jidprimer,
                        inside, ref_const, universal_rev, reverse_primer_group,
                        species)
    check_barcodes(primer_tag_file_count, primer_tag_file, Fail_file, sample,
                   Output_trim, threshold_barcode)
    print_trimmed_sequences(Output_trim, primer_tag_file,
                            primer_tag_file_count, inside, ref_const)
    return ()


def dfgdfg():
    out, ind = '', 0
    for s in seqs:
        ind = 0
        if (len(seqs[s]) > 1):
            for id in seqs[s]:
                ind = ind + 1
                if (ind == 1):
                    if (id.count("|") != 0):
                        f = [0] * len(id.split("|")[1].split("_"))
                    else:
                        f = 0
                freq = list(
                    map(int,
                        id.split("__")[1].split("|")[0].split("_")))
                if (len(freq) > 1):
                    f = list(map(np.add, f, freq))
                else:
                    f = f + freq[0]
            if (id.count("|") != 0):
                id = id.split("__")[0] + "__" + "_".join(map(
                    str, f)) + "|" + id.split("|")[1]
            else:
                id = id.split("__")[0] + "__" + str(f)
            out = out + ">" + id + "\n" + s + "\n"
        else:
            for id in seqs[s]:
                break
            out = out + ">" + id + "\n" + s + "\n"
        ind = ind + 1
        if (ind > 500):
            write_out(out, nn_orf_filtered)
            out, ind = '', 0
    write_out(out, nn_orf_filtered)
    out, ind = '', 0
    print("\nTotal sequences:", number_seqs,
          "\tNumber passed with open reading frame:", number_accepted,
          "\tPercentage accepted:", number_accepted * 100.0 / number_seqs, "%")
    return ()


def check_barcodes(primer_tag_file_count, primer_tag_file, Fail_file, sample,
                   Output_trim, threshold):
    # mode = "CONTINUE"
    mode = "INIT"
    if (mode == "CONTINUE"):
        bc_complete = get_barcodes_completed(primer_tag_file)
    elif (mode == "INIT"):
        bc_complete = {}
        outs = [
            '',
            ('#ID\tnumber_of_reads\ttotal_reads_with_BC\tJ_barcode\tV_barcode\t',
             'bp_mismatches_to_consensus\tBC_accepted\tconsensus\tsequence\n')
        ]
        files = [Output_trim, primer_tag_file]
        for i in range(0, len(files)):
            fh = open(files[i], "w")
            fh.write(outs[i])
            fh.close()
    fh = open(primer_tag_file_count, "r")
    seqs = Tree()
    for l in fh:
        if (l[0] != "#"):
            l = l.strip().split()
            v_tag, j_tag, sequence, header = l[2], l[1], l[3], l[0]
            if (j_tag + "\t" + v_tag not in bc_complete):
                seqs[j_tag + "\t" + v_tag][sequence][header].value = 1
    fh.close()
    print(len(seqs), "Unique tags")
    failed, ind_f, out, outp, ind = '', 0, '', '', 0
    unique_seqs = Tree()
    min_depth = (1.0 / (1 - threshold))
    # indp, total_tags, passed_seqs_total = 0, 0, 0
    total_tags, passed_seqs_total = 0, 0
    thresh = 20
    detailed = "FALSE"  # if detailed output required
    for t1 in seqs:
        # t, u_seq, u_freq, u_header = 0, [], [], []
        u_seq, u_freq, u_header = [], [], []
        for s in seqs[t1]:
            f = 0
            for h in seqs[t1][s]:
                f = f + int(h.split(":")[1])
                break
            total_tags = total_tags + f
            u_seq.append(s)
            u_freq.append(f)
            u_header.append(h)
            ind = ind + 1
            if (len(u_freq) > 500):
                break
        f = sum(u_freq)
        h = h.split(":")[0].split("__")[0] + "__" + str(sum(u_freq))
        if (len(u_freq) > 100):
            print("BC group size:\t", len(u_freq))
        if (len(seqs[t1]) == 1):
            passed_seqs_total = passed_seqs_total + f
            unique_seqs[s][h][t1].value = 1
            outp = outp + h + "\t" + str(f) + "\t" + str(
                f) + "\t" + t1 + "\t0\tYES\t" + s + "\t" + s + "\n"
            if (f > thresh):
                print(t1, f)
        else:
            if (
                    sum(u_freq) >= min_depth
            ):  # if fewer sequences associated with J-barcode, but not all are identical, then reject
                if (max(u_freq) * 1.0 / sum(u_freq) > threshold):
                    print("\t\t1:", t1, sum(u_freq), min_depth)
                    s = u_seq[u_freq.index(max(u_freq))]
                    unique_seqs[s][h][t1].value = 1
                    passed_seqs_total = passed_seqs_total + f
                    if (detailed == "TRUE"):
                        mm = []
                        l1 = len(consensus)
                        for i in range(0, len(u_seq)):
                            s1, s2, p1 = trim_sequences(
                                consensus, u_seq[i], l1, len(u_seq[i]))
                            if (p1 == 1):
                                p, mismatches = get_diff(s1, s2, 4)
                            else:
                                mismatches = 5
                            mm.append(mismatches)
                            mismatches = get_mismatches_from_consensus(
                                u_seq[i], s)
                            # if(mismatches>4):mismatches=">4"
                            outp = (outp + u_header[i] + "\t" +
                                    str(u_freq[i]) + "\t" + str(sum(u_freq)) +
                                    "\t" + t1 + "\t" + str(mismatches) +
                                    "\tYES\t" + s + "\t" + u_seq[i] + "\n")
                            if (i > 200):
                                break
                    else:
                        outp = outp + h + "\tNA\t" + str(sum(
                            u_freq)) + "\t" + t1 + "\tNA\tYES\t" + s + "\tNA\n"
                    if (f > thresh):
                        print(t1, f)
                else:
                    print("\t\t2:", t1, sum(u_freq), min_depth)
                    print(tmp_file)
                    consensus = get_consensus_sequence(u_seq, u_freq, tmp_file,
                                                       threshold)
                    if (consensus.count("_") == 0 and len(consensus) > 3):
                        unique_seqs[consensus][h][t1].value = 1
                        passed_seqs_total = passed_seqs_total + f
                        # if(detailed == "TRUE"):
                        mm = []
                        l1 = len(consensus)
                        for i in range(0, len(u_seq)):
                            s1, s2, p1 = trim_sequences(
                                consensus, u_seq[i], l1, len(u_seq[i]))
                            if (p1 == 1):
                                p, mismatches = get_diff(s1, s2, 4)
                            else:
                                mismatches = 5
                            # mismatches = get_mismatches_from_consensus(u_seq[i],consensus)
                            mm.append(mismatches)
                        if (mm.count(0) >= 2):
                            for i in range(0, len(u_seq)):
                                outp = (outp + u_header[i] + "\t" +
                                        str(u_freq[i]) + "\t" +
                                        str(sum(u_freq)) + "\t" + t1 + "\t" +
                                        str(mm[i]) + "\tYES\t" + consensus +
                                        "\t" + u_seq[i] + "\n")
                            else:
                                outp = (outp + u_header[i] + "\t" +
                                        str(u_freq[i]) + "\t" +
                                        str(sum(u_freq)) + "\t" + t1 + "\t" +
                                        str(mm[i]) + "\tNO\t" + consensus +
                                        "\t" + u_seq[i] + "\n")
                    else:
                        failed, ind_f = failed + ">" + h + "_BC_multiplicity\n" + s + "\n", ind_f + 1
                        if (detailed == "TRUE"):
                            for i in range(0, len(u_seq)):
                                print(t1, "\t", i)
                                mismatches = "NA"
                                outp = outp + u_header[i] + "\t" + str(
                                    u_freq[i]
                                ) + "\t" + str(
                                    sum(u_freq)
                                ) + "\t" + t1 + "\t" + mismatches + "\tNO\t" + consensus + "\t" + u_seq[
                                    i] + "\n"
            else:
                consensus = get_consensus_sequence(u_seq, u_freq, tmp_file,
                                                   threshold)
                if (consensus.count("_") == 0 and len(consensus) > 3):
                    mm = []
                    l1 = len(consensus)
                    for i in range(0, len(u_seq)):
                        s1, s2, p1 = trim_sequences(consensus, u_seq[i], l1,
                                                    len(u_seq[i]))
                        if (p1 == 1):
                            p, mismatches = get_diff(s1, s2, 4)
                        else:
                            mismatches = 5
                        # mismatches = get_mismatches_from_consensus(u_seq[i],consensus)
                        mm.append(mismatches)
                    # print "less than 5 mm",mm
                    if (mm.count(0) >= 2):
                        unique_seqs[consensus][h][t1].value = 1
                        passed_seqs_total = passed_seqs_total + f
                        for i in range(0, len(u_seq)):
                            outp = outp + h + "\tNA\t" + str(
                                sum(u_freq[i])
                            ) + "\t" + t1 + "\tNA\tYES\t" + consensus + "\t" + u_seq[
                                i] + "\n"
                        # if(mismatches>4):mismatches=">4"
                        # outp=outp+u_header[i]+"\t"+str(u_freq[i])+"\t"+str(sum(u_freq))+"\t"+t1+"\t"+str(mismatches)+"\tYES\t"+consensus+"\t"+u_seq[i]+"\n"
                    # else:outp=outp+h+"\tNA\t"+str(sum(u_freq))+"\t"+t1+"\tNA\tYES\t"+consensus+"\tNA\n"
                    # if(f>thresh):print t1, f
                    else:
                        failed, ind_f = failed + ">" + h + "_BC_multiplicity\n" + consensus + "\n", ind_f + 1
                        print("FAILED")
                else:
                    failed, ind_f = failed + ">" + h + "_BC_multiplicity\n" + consensus + "\n", ind_f + 1
                    print(+h + "_BC_multiplicity\n" + consensus)
                    for i in range(0, len(u_seq)):
                        mismatches = "NA"
                        outp = outp + u_header[i] + "\t" + str(
                            u_freq[i]
                        ) + "\t" + str(
                            sum(u_freq)
                        ) + "\t" + t1 + "\t" + mismatches + "\tNO\t" + consensus + "\t" + u_seq[
                            i] + "\n"
        if (ind > 100):
            write_out(outp, primer_tag_file)
            outp, ind = '', 0
        if (ind_f > 500):
            write_out(failed, Fail_file)
            failed, ind_f = '', 0
    write_out(outp, primer_tag_file)
    write_out(failed, Fail_file)
    print("Getting refined sequences")
    del failed, ind_f, out, seqs, outp
    return ()


def dfg():
    for i in h:
        print(len(unique_seqs))
        out, ind = '', 0
        reg = []
        for i in regions:
            print(i)
            reg.append(i)
        reg.sort()
        header = "_".join(reg)
        for seq in unique_seqs:
            f = [0] * len(reg)
            for const in unique_seqs[seq]:
                ind1 = reg.index(const)
                f[ind1] = f[ind1] + len(unique_seqs[seq][const])
            print(f)
            for id in unique_seqs[seq][const]:
                break
            out = out + ">" + id.split("__")[0] + "__" + "_".join(map(
                str, f)) + "|" + header + "\n" + seq + "\n"
            ind = ind + 1
            if (ind > 100):
                write_out(out, Output_trim)
                out, ind = '', 0
        write_out(out, Output_trim)
    else:
        print("EDIT CODE HERE: print out trimmed sequences")
