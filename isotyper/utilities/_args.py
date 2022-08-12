#!/usr/bin/env python
from pathlib import Path
import sys

# from arguments
SAMPLE_ID = sys.argv[1]  # name of sample. Also the prefix of the file.
SOURCE = Path(sys.argv[2])
OUT_PATH = Path(sys.argv[3])
ORG = sys.argv[4]  # organism
LENGTH = sys.argv[5]
COMMAND_MODE = sys.argv[6]
# output subfolders
OUT_FASTQ = OUT_PATH / "FASTQ_FILES"
OUT_ORTSEQ = OUT_PATH / "ORIENTATED_SEQUENCES"
OUT_ORTSEQ_TMP = OUT_ORTSEQ / "TMP"
OUT_NET = OUT_ORTSEQ / "NETWORKS"
