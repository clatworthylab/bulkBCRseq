#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from igraph import *


def main(argv):
    """Summary

    Parameters
    ----------
    argv : TYPE
        Description
    """
    if len(argv) < 2:
        print("usage: python %s sample directory" % argv[0])
        sys.exit(1)

    sample = argv[1]

    if len(argv) == 2:
        directory = "./"
    else:
        directory = argv[2]
        if not directory.endswith("/"):
            directory = directory + "/"

    att_file = directory + "Att_" + sample + ".txt"
    out_file = directory + "Formatted_Att_" + sample + ".txt"
    cluster_file = directory + "Plot_ids_" + sample + ".txt"

    att = pd.read_csv(att_file, sep="\t", header=None)
    cluster = pd.read_csv(cluster_file, header=None, skiprows=1, sep="\t")
    cluster_id = cluster[0]
    cluster_value = cluster_id.value_counts()
    cluster_dict = {}
    i = 1
    for c, v in cluster_value.items():
        cluster_dict[c] = i
        i += 1
    cluster[4] = [cluster_dict[c] for c in cluster[0]]
    cluster_id_dict = dict(zip(cluster[1], cluster[0]))
    cluster_sorted_dict = dict(zip(cluster[1], cluster[4]))

    att[3] = [cluster_id_dict[c] for c in att[0]]

    att[4] = [cluster_sorted_dict[c] for c in att[0]]

    top_10_clusters = [str(x) for x in list(cluster_value.index[:10])]
    top_10_clusters_dict = dict(
        zip(top_10_clusters, ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"])
    )
    att[5] = [str(x) for x in cluster[3]]
    for x in att.index:
        if att.at[x, 5] not in top_10_clusters:
            att.at[x, 5] = "K"
        else:
            att.at[x, 5] = top_10_clusters_dict[att.at[x, 5]]

    isotype = []
    for c in att[0]:
        isotype.append(
            [int(x) for x in c.split("__")[1].split("|")[0].split("_")]
        )
    isotype = pd.DataFrame(
        isotype, columns=[c.split("__")[1].split("|")[1].split("_")]
    )
    isotypes = list(isotype.idxmax(axis=1))
    isotypes = [str(i).strip("()|,|''") for i in isotypes]
    isotype_dict = dict(zip(att[0], isotypes))
    att[6] = [isotype_dict[c] for c in att[0]]

    if "-_-" in att[0][0]:
        tissue = []
        tissue_conversion = {
            "bm": "bone_marrow",
            "l": "lung",
            "t": "tongue",
            "b": "bladder",
            "m": "meninges",
            "ln": "mesenteric_lymph_nodes",
            "d": "duodenum",
            "i": "ilieum",
            "cp": "colon_proximal",
            "cd": "colon_distal",
        }

        for c in att[0]:
            tissue.append(c.split("-_-")[0].split("_M")[1][1:])
        tissue_dict = dict(zip(att[0], tissue))
        att[7] = [tissue_conversion[tissue_dict[c]] for c in att[0]]

    att.to_csv(out_file, sep="\t", header=None, index=False)

    node_file = directory + "Formatted_Att_" + sample + ".txt"
    edge_file = directory + "Edges_" + sample + ".txt"

    nodes = pd.read_csv(node_file, sep="\t", header=None)
    edge = pd.read_csv(edge_file, sep="\t", header=None)

    edges = [tuple(e) for e in edge.values]
    edges = [e for e in edges if e[0] > e[1]]
    edges_list = [tuple((e[0], e[1])) for e in edges]
    g = Graph.Formula()
    g.add_vertices([n for n in nodes[0]])
    g.add_edges(edges_list)
    g.vs["size"] = [int(s) for s in nodes[1]]
    remove_ids1 = [v.index for v in g.vs if v.degree() == 0]
    remove_ids2 = [v.index for v in g.vs if v["size"] == 1]
    remove_ids = list(set(remove_ids1 + remove_ids2))

    clone_dict = dict(zip(nodes[0], nodes[5]))
    iso_dict = dict(zip(nodes[0], nodes[6]))
    e_clone = [clone_dict[e[0]] for e in edges]
    e_iso = [iso_dict[e[0]] for e in edges]
    g.vs["clone"] = nodes[5]
    g.vs["isotype"] = nodes[6]

    if nodes.shape[1] > 7:
        tis_dict = dict(zip(nodes[0], nodes[7]))
        e_tis = [tis_dict[e[0]] for e in edges]
        g.vs["tissue"] = nodes[7]

    if nodes.shape[1] > 7:
        e_anno = {0: e_clone, 1: e_iso, 2: e_tis}
        e_df = pd.DataFrame(e_anno)
        g.es["clone"] = e_df[0]
        g.es["isotype"] = e_df[1]
        g.es["tissue"] = e_df[2]
    else:
        e_anno = {0: e_clone, 1: e_iso}
        e_df = pd.DataFrame(e_anno)
        g.es["clone"] = e_df[0]
        g.es["isotype"] = e_df[1]

    g.delete_vertices(remove_ids)
    max_freq = max([int(v["size"]) for v in g.vs])

    lyt = g.layout_graphopt(niter=800)

    visual_style = {}
    visual_style["vertex_size"] = [1 + s * 10 / max_freq for s in g.vs["size"]]
    visual_style["vertex_frame_width"] = 0.1
    visual_style["vertex_label"] = g.vs["name"]
    visual_style["vertex_label_size"] = 0
    visual_style["layout"] = lyt
    visual_style["bbox"] = (800, 800)
    visual_style["margin"] = 20
    visual_style["inline"] = True

    isotype_col_dict = {
        "IGHA": "#4e79a7",
        "IGHD": "#f28e2b",
        "IGHE": "#e15759",
        "IGHG1": "#76b7b2",
        "IGHG2A": "#59a14f",
        "IGHG2B": "#edc948",
        "IGHG2C": "#b07aa1",
        "IGHG3": "#ff9da7",
        "IGHM": "#9c755f",
        "NaN": "#f2f2f2",
        np.nan: "#f2f2f2",
    }
    tissue_col_dict = {
        "bone_marrow": "#4E79A7",
        "lung": "#A0CBE8",
        "tongue": "#F28E2B",
        "bladder": "#FFBE7D",
        "meninges": "#59A14F",
        "mesenteric_lymph_nodes": "#8CD17D",
        "duodenum": "#B6992D",
        "ilieum": "#F1CE63",
        "colon_proximal": "#499894",
        "colon_distal": "#86BCB6",
    }
    clone_col_dict = {
        "A": "#4e79a7",
        "B": "#a0cbe8",
        "C": "#f28e2b",
        "D": "#ffbe7d",
        "E": "#59a14f",
        "F": "#8cd17d",
        "G": "#b6992d",
        "H": "#f1ce63",
        "I": "#e15759",
        "J": "#ff9d9a",
        "K": "#f2f2f2",
    }

    # file names
    outI = directory + sample + "_isotype.png"
    outC = directory + sample + "_top_10_clones.png"
    outT = directory + sample + "_tissue.png"

    visual_style["edge_color"] = [isotype_col_dict[i] for i in g.es["isotype"]]
    visual_style["vertex_color"] = [
        isotype_col_dict[i] for i in g.vs["isotype"]
    ]
    visual_style["vertex_frame_color"] = [
        isotype_col_dict[i] for i in g.vs["isotype"]
    ]
    plot(g, outI, **visual_style)

    visual_style["edge_color"] = [clone_col_dict[i] for i in g.es["clone"]]
    visual_style["vertex_color"] = [clone_col_dict[i] for i in g.vs["clone"]]
    visual_style["vertex_frame_color"] = [
        clone_col_dict[i] for i in g.vs["clone"]
    ]
    plot(g, outC, **visual_style)

    if nodes.shape[1] > 7:
        visual_style["edge_color"] = [
            tissue_col_dict[i] for i in g.es["tissue"]
        ]
        visual_style["vertex_color"] = [
            tissue_col_dict[i] for i in g.vs["tissue"]
        ]
        visual_style["vertex_frame_color"] = [
            tissue_col_dict[i] for i in g.vs["tissue"]
        ]
        plot(g, outT, **visual_style)


if __name__ == "__main__":
    main(sys.argv)
