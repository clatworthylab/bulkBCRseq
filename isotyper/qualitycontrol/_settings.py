#!/usr/bin/env python
from isotyper.utilities._settings import LIBPATH
from isotyper.utilities._args import (
    ORG,
    OUT_FASTQ,
    OUT_NET,
    OUT_ORTSEQ,
    OUT_ORTSEQ_TMP,
    OUT_PATH,  # for import into main script
    SAMPLE_ID,
)

# primer file
PRIMER_FILE = LIBPATH / f"Primers_{ORG}_IG_RBR_Constant_region_MPLX_table.txt"
CODON_FILE = LIBPATH / f"Codon_table2.txt"

# reference files
REFJ = LIBPATH / f"Reference_nn_{ORG}_IGHJ.fasta"
REFVP = LIBPATH / f"Reference_protein_{ORG}_IGHV.fasta"
REF_CONST = LIBPATH / f"Reference_nn_{ORG}_IGH_constant_exon1.fasta"

# get overlapping paired reads
SEQ_FASTA_FILE1 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_1.fasta"
SEQ_FASTA_FILE2 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_2.fasta"
UNTRIMMED_FASTA = OUT_ORTSEQ_TMP / f"Untrimmed_{SAMPLE_ID}.fasta"

# trimming
TMP_FILE = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.txt"
FAIL_FILE = OUT_FASTQ / f"Fail_filtered_{SAMPLE_ID}.fasta"
PRIMER_TAG_FILE = (
    OUT_ORTSEQ_TMP / f"Barcode_filtering_information_{SAMPLE_ID}.txt"
)
PRIMER_TAG_FILE_COUNT = OUT_ORTSEQ_TMP / f"All_barcodes_{SAMPLE_ID}.txt"

# trimmed
TRIM1_ALL = OUT_ORTSEQ_TMP / f"trimmed_orientated_all_{SAMPLE_ID}.fasta"
TRIM2_J = OUT_ORTSEQ_TMP / f"Filtered_J_{SAMPLE_ID}.fasta"
TRIM3_RED = OUT_ORTSEQ_TMP / f"Filtered_reduced_{SAMPLE_ID}.fasta"

# filtered
FILTERED_OUT = OUT_ORTSEQ / f"Filtered_ORFs_sequences_all_{SAMPLE_ID}.fasta"
FILTERED_OUT_NT = OUT_ORTSEQ / f"Nucleotide_ORF_filtered_all_{SAMPLE_ID}.fasta"
FILTERING_REPORT = OUT_ORTSEQ / f"Filtering_report_{SAMPLE_ID}.txt"
