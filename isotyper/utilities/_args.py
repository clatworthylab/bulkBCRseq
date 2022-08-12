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
NCPUS = sys.argv[7]
# output subfolders
OUT_FASTQ = OUT_PATH / "FASTQ_FILES"
OUT_ORTSEQ = OUT_PATH / "ORIENTATED_SEQUENCES"
OUT_ORTSEQ_TMP = OUT_ORTSEQ / "TMP"
OUT_NET = OUT_ORTSEQ / "NETWORKS"
# FRSEQ_FILE = OUT_NET / f"Fully_reduced_{SAMPLE_ID}.fasta
CLUST_ID_FILE = OUT_NET / f"Cluster_identities_{SAMPLE_ID}.txt"
# statistic results
NETSTATS = OUT_PATH / f"Network_statistics_{SAMPLE_ID}.txt"
CLUSTER_SIZE_DIST = OUT_PATH / f"Distribution_cluster_sizes_{SAMPLE_ID}.txt"
VERTEX_SIZE_DIST = OUT_PATH / f"Distribution_vertex_sizes_{SAMPLE_ID}.txt"
NETSTATS_PER_CHAIN = OUT_PATH / f"Network_statistics_per_chain_{SAMPLE_ID}.txt"
