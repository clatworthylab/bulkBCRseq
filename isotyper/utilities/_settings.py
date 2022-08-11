#!/usr/bin/env python
from pathlib import Path
from isotyper import __path__
from isotyper.utilities._args import ORG, OUT_NET, SAMPLE_ID

ISOTYPERPREFIX = Path(__path__[0])
# set up some default paths
LIBPATH = ISOTYPERPREFIX / "library"
EXTPATH = ISOTYPERPREFIX / "external"
PERLCMD = "$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}"

EDGE_LENGTHS = 0.85
READ_NUMBER_DIVISION = "__"
THRESHOLD_BARCODE = 0.80  # Certainty for accepting a barcode
THRESHOLD_GENE_SCORE = 5 / 2
REVERSE_COMPLEMENT_DICT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    ".": ".",
}
MIN_QUAL = 32
# change here if necessary
R1PATTERN = "_R1_001"
R2PATTERN = "_R2_001"
PRIMER_FILE = LIBPATH / f"Primers_{ORG}_IG_RBR_Constant_region_MPLX_table.txt"
TMP_FILE = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.txt"
