#!/usr/bin/env python

from isotyper.utilities._args import OUT_NET, OUT_PATH, SAMPLE_ID

OUT_STAT = OUT_PATH / "STATISTICS"
# FRSEQ_FILE = OUT_NET / f"Fully_reduced_{SAMPLE_ID}.fasta
CLUST_ID_FILE = OUT_NET / f"Cluster_identities_{SAMPLE_ID}.txt"
# statistic results
NETSTATS = OUT_STAT / f"Network_statistics_{SAMPLE_ID}.txt"
CLUSTER_SIZE_DIST = OUT_STAT / f"Distribution_cluster_sizes_{SAMPLE_ID}.txt"
VERTEX_SIZE_DIST = OUT_STAT / f"Distribution_vertex_sizes_{SAMPLE_ID}.txt"
NETSTATS_PER_CHAIN = OUT_STAT / f"Network_statistics_per_chain_{SAMPLE_ID}.txt"
