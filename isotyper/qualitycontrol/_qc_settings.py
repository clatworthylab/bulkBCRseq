#!/usr/bin/env python
from isotyper.utilities._settings import LIBPATH
from isotyper.utilities._args import (
    ORG,
    SAMPLE_ID,
    OUT_NET,
    OUT_FASTQ,
    OUT_ORTSEQ_TMP,
)


# Files for Step2
PRIMER_FILE = LIBPATH / f"Primers_{ORG}_IG_RBR_Constant_region_MPLX_table.txt"
TMP_FILE = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.txt"
FAIL_FILE = OUT_FASTQ / f"Fail_filtered_{SAMPLE_ID}.fasta"
SEQ_FASTA_FILE1 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_1.fasta"
SEQ_FASTA_FILE2 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_2.fasta"
UNTRIMMED_FASTA = OUT_ORTSEQ_TMP / f"Untrimmed_{SAMPLE_ID}.fasta"
