#!/usr/bin/env python
from pathlib import Path
import sys

# from arguments
OUT_PATH = Path(sys.argv[1])
SAMPLE_ID = sys.argv[2]  # name of sample. Also the prefix of the file.
ORG = sys.argv[6]  # organism
SOURCE = Path(sys.argv[7])
LENGTH = sys.argv[8]
COMMAND_MODE = sys.argv[11]
# if len(sys.argv) > 13:
#     REVERSE_PRIMER_GROUP = sys.argv[13]
# else:
#     REVERSE_PRIMER_GROUP = "OTHER"
# output subfolders
OUT_FASTQ = OUT_PATH / "FASTQ_FILES"
OUT_ORTSEQ = OUT_PATH / "ORIENTATED_SEQUENCES"
OUT_ORTSEQ_TMP = OUT_ORTSEQ / "TMP"
OUT_NET = OUT_ORTSEQ / "NETWORKS"
