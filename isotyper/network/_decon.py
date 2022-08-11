#!/usr/bin/env python
from operator import add
from pathlib import Path
from typing import Dict, Tuple, List

from isotyper.utilities import (
    fasta_iterator,
    write_out,
    Tree,
    deconvolute_same_array,
    READ_NUMBER_DIVISION,
)
from isotyper.network._settings import CLUSTERED_TMP_FILE


def deconvolute_edges(
    seq_file: Path,
    ver_rel_file: Path,
    att_file: Path,
    file_seqs: Path,
    file_edges: Path,
):
    """Deconvolute edges.

    Parameters
    ----------
    seq_file : Path
        path to filtered nucleotide sequences.
    ver_rel_file : Path
        path to vertex relation file.
    att_file : Path
        path to att file.
    file_seqs : Path
        path to sequences file.
    file_edges : Path
        path to edges file
    """
    decon_identical(
        seq_file=seq_file,
        ver_rel_file=ver_rel_file,
        att_file=att_file,
        file_seqs=file_seqs,
    )
    decon_edges(
        ver_rel_file=ver_rel_file,
        file_seqs=file_seqs,
        file_edges=file_edges,
        att_file=att_file,
        seq_file=seq_file,
    )
    reduce_edges(file_in=CLUSTERED_TMP_FILE, file_out=file_edges)


def print_vertices(
    all_seqs: Dict,
    inverse: Dict,
    same: Dict,
    att_file: Path,
    length: Dict,
    info: str,
    freq_id: Dict,
):
    """Summary

    Parameters
    ----------
    all_seqs : Dict
        Description
    inverse : Dict
        Description
    same : Dict
        Description
    att_file : Path
        Description
    length : Dict
        Description
    info : str
        Description
    freq_id : Dict
        Description
    """
    out = ""
    # total = 0
    fh = open(att_file, "w")
    fh.close()
    ind = 0
    for seq_id in all_seqs:
        ind = ind + 1
        if seq_id in inverse:
            if seq_id in same:
                id1 = (
                    seq_id.split(READ_NUMBER_DIVISION)[0]
                    + READ_NUMBER_DIVISION
                    + "_".join(map(str, same[seq_id]))
                    + info
                )
                if seq_id not in length:
                    out = (
                        out
                        + id1
                        + "\t"
                        + str(same[seq_id])
                        + "\t"
                        + all_seqs[seq_id]
                        + "\n"
                    )
                else:
                    out = (
                        out
                        + id1
                        + "\t"
                        + str(sum(same[seq_id]))
                        + "\t"
                        + length[seq_id]
                        + "\n"
                    )
        else:
            if seq_id in freq_id:
                id1 = freq_id[seq_id]
            else:
                print(seq_id)
            freq = sum(
                map(
                    int,
                    id1.split(READ_NUMBER_DIVISION)[1].split("|")[0].split("_"),
                )
            )
            out = out + id1 + "\t" + str(freq) + "\t" + all_seqs[seq_id] + "\n"
        if ind > 300:
            write_out(out, att_file)
            ind = 0
            out = ""
    write_out(out, att_file)


def decon_identical(
    seq_file: Path,
    ver_rel_file: Path,
    att_file: Path,
    file_seqs: Path,
):
    """Summary

    Parameters
    ----------
    seq_file : Path
        Description
    ver_rel_file : Path
        Description
    att_file : Path
        Description
    file_seqs : Path
        Description
    """
    fh = open(seq_file, "r")
    all_seqs, seqs, freq_id = {}, {}, {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split(READ_NUMBER_DIVISION)[0]] = sequence
        all_seqs[header.split(READ_NUMBER_DIVISION)[0]] = sequence
        freq_id[header.split(READ_NUMBER_DIVISION)[0]] = header
    fh.close()
    fh1 = open(ver_rel_file, "r")
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
                        j1.split(READ_NUMBER_DIVISION)[1]
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
            .split(READ_NUMBER_DIVISION)[1]
            .split("|")
        )
        >= 2
    ):
        info = (
            "|"
            + freq_id[j.split(READ_NUMBER_DIVISION)[0]]
            .split(READ_NUMBER_DIVISION)[1]
            .split("|")[1]
        )
    write_out(out, file_seqs)
    print_vertices(
        all_seqs=all_seqs,
        inverse=inverse,
        same=same,
        att_file=att_file,
        length=length,
        info=info,
        freq_id=freq_id,
    )


def get_inverse_ids(
    file_seqs: Path, att_file: Path, freq_id: Dict
) -> Tuple[Dict, Dict]:
    """Summary

    Parameters
    ----------
    file_seqs : Path
        Description
    att_file : Path
        Description
    freq_id : Dict
        Description

    Returns
    -------
    Tuple[Dict, Dict]
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
    fh = open(att_file, "r")
    for l in fh:
        l = l.strip().split()
        if l[0].split(READ_NUMBER_DIVISION)[0] in raw:
            raw[l[0].split(READ_NUMBER_DIVISION)[0]] = l[0]
    fh.close()
    return (inverse, raw)


def print_single_edges(
    file_edges: Path, inverse: Dict, edges: Tree, tmp_file: Path, raw: Dict
):
    """Summary

    Parameters
    ----------
    file_edges : Path
        Description
    inverse : Dict
        Description
    edges : Tree
        Description
    tmp_file : Path
        Description
    raw : Dict
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


def decon_edges(
    ver_rel_file: Path,
    file_seqs: Path,
    file_edges: Path,
    att_file: Path,
    seq_file: Path,
):
    """Summary

    Parameters
    ----------
    ver_rel_file : Path
        Description
    file_seqs : Path
        Description
    file_edges : Path
        Description
    att_file : Path
        Description
    seq_file : Path
        Description
    """
    fh = open(seq_file, "r")
    freq_id = {}
    for header, sequence in fasta_iterator(fh):
        freq_id[header.split(READ_NUMBER_DIVISION)[0]] = header
    fh.close()
    inverse, raw = get_inverse_ids(
        file_seqs=file_seqs, att_file=att_file, freq_id=freq_id
    )
    edges = Tree()
    fh1 = open(ver_rel_file, "r")
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
    print_single_edges(
        file_edges=file_edges,
        inverse=inverse,
        edges=edges,
        tmp_file=CLUSTERED_TMP_FILE,
        raw=raw,
    )


def reduce_edges(file_in: Path, file_out: Path):
    """Reduce edges.

    Parameters
    ----------
    file_in : Path
        Description
    file_out : Path
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
