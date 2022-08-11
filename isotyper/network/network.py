#!/usr/bin/env python
from pathlib import Path
from typing import Dict

from isotyper.utilities import fasta_iterator
from isotyper.utilities import READ_NUMBER_DIVISION


def get_seqs_single(file: Path) -> Dict:
    """return dictionary of sequences from fasta file.

    Parameters
    ----------
    file : Path
        fasta file.

    Returns
    -------
    Dict
        sequences
    """
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[header.split(READ_NUMBER_DIVISION)[0]] = sequence
    fh.close()
    return seqs
