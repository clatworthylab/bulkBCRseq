#!/usr/bin/env python
import os
import pandas as pd
import sys


def fasta_iterator(fh):
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
    while True:
        line = fh.readline()
        if line.startswith(">"):
            break
    while True:
        header = line[1:-1].rstrip()
        sequence = fh.readline().rstrip()
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith(">"):
                break
            sequence += line.rstrip()
        yield (header, sequence)
        if not line:
            return


def Write_output(out, file):
    """Summary

    Parameters
    ----------
    out : TYPE
        Description
    file : TYPE
        Description
    """
    fh = open(file, "a")
    fh.write(out)
    fh.close()


def main(argv):
    """Summary

    Parameters
    ----------
    argv : TYPE
        Description
    """
    if len(argv) < 2:
        print("usage: python %s sample directory downsample" % argv[0])
        sys.exit(1)

    fasta_path = argv[1]
    airr_path = argv[2]
    if len(argv) > 3:
        downsample = argv[3]

    fh = open(fasta_path, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[header] = sequence
    fh.close()

    df = pd.read_csv(airr_path, sep="\t")
    if "downsample" in locals():
        df = df.sample(int(downsample))
        df.to_csv(airr_path.replace(".tsv", "_downsample.tsv"), sep="\t")
    df["sequence_id"] = [
        s + "-_-" + f for s, f in zip(df["sample_id"], df["sequence_id"])
    ]
    headers = []
    headers = [h for h in df["sequence_id"]]

    seqs2 = {
        filtered_header: seqs[filtered_header] for filtered_header in headers
    }
    if "downsample" in locals():
        out_file = airr_path.replace(".tsv", "_downsample_fully_reduced.fasta")
    else:
        out_file = airr_path.replace(".tsv", "_fully_reduced.fasta")

    fh1 = open(out_file, "w")
    fh1.close()
    out = ""
    for l in seqs2:
        out = ">" + l + "\n" + seqs2[l] + "\n"
        Write_output(out, out_file)


if __name__ == "__main__":
    main(sys.argv)
