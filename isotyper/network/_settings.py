#!/usr/bin/env python
from isotyper.utilities._args import (
    SAMPLE_ID,
    OUT_NET,
)

# clustered files
CLUSTERED_TMP_FILE = OUT_NET / f"NN_Tmp_cluster_{SAMPLE_ID}.1"
COCLUSTERED = CLUSTERED_TMP_FILE.with_suffix("._coclustered")
# network files
VERTEX_REL_FILE = OUT_NET / f"Vertex_relations_{SAMPLE_ID}.txt"
ATT_FILE = OUT_NET / f"Att_{SAMPLE_ID}.txt"
SEQS_FILE = OUT_NET / f"Sequences_{SAMPLE_ID}.txt"
EDGES_FILE = OUT_NET / f"Edges_{SAMPLE_ID}.txt"
CLUSTER_FILE = OUT_NET / f"Cluster_identities_{SAMPLE_ID}.txt"
# final output
REDUCED_FILE = OUT_NET / f"Fully_reduced_{SAMPLE_ID}.fasta"
# other output
PLOT_IDS_FILE = OUT_NET / f"Plot_ids_{SAMPLE_ID}.txt"
CHKEDGES_FILE = OUT_NET / f"Checked_edges_{SAMPLE_ID}.txt"
