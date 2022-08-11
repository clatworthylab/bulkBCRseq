#!/usr/bin/env python
from isotyper.utilities._args import (
    ORG,
    SAMPLE_ID,
    OUT_NET,
    OUT_FASTQ,
    OUT_ORTSEQ,
    OUT_ORTSEQ_TMP,
)

# clustered file
# tmp_file_1 = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.1"
CLUSTERED_TMP_FILE = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.1"
