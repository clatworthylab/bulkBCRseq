#!/usr/bin/env python
from operator import itemgetter
from typing import Tuple, List

from isotyper.utilities._args import *
from isotyper.utilities import write_out, create_file, Tree


def renyi_entropy(
    cpoints: List,
    cvdf: List,
    vpoints: List,
    vvdf: List,
    totalv: int,
    totalreads: int,
) -> Tuple[float, float]:
    """Renyi entropy

    Parameters
    ----------
    cpoints : List
        cluster sizes
    cvdf : List
        cluster size distribution
    vpoints : List
        vertex sizes
    vvdf : List
        vertex size distribution
    totalv : int
        total vertices
    totalreads : int
        total read counts

    Returns
    -------
    Tuple[float, float]
        vertex and cluster Renyi entropies.
    """
    vrenyi = 0
    crenyi = 0
    tv = totalreads * totalreads * 1.0  # KT: what is this * 1.0 for?
    tc = totalv * totalv * 1.0
    for i in range(0, len(vpoints)):
        vrenyi += vvdf[i] * (vpoints[i] * vpoints[i] / tv)
    for i in range(0, len(cpoints)):
        crenyi += cvdf[i] * (cpoints[i] * cpoints[i] / tc)
    return (vrenyi, crenyi)


def vdf(n: List) -> Tuple[List, List]:
    """distribution function.

    Parameters
    ----------
    n : List
        list of numbers

    Returns
    -------
    Tuple[List, List]
        size vector and cumulative distribution function
    """
    points = sorted(list(set(n)))
    vdf = []
    for i in range(0, len(points)):
        vdf.append(n.count(points[i]))
    return (points, vdf)


def gini_index(
    cpoints: List, cvdf: List, vpoints: List, vvdf: List
) -> Tuple[float, float]:
    """Gini index caculator for vertex and cluster size.

    Parameters
    ----------
    cpoints : List
        cluster sizes
    cvdf : List
        cluster size distribution
    vpoints : List
        vertex sizes
    vvdf : List
        vertex size distribution

    Returns
    -------
    Tuple[float, float]
        vertex and cluster gini indices.
    """
    (vgini) = get_gini(n=vpoints, v=vvdf)
    (cgini) = get_gini(n=cpoints, v=cvdf)
    return (vgini, cgini)


def get_gini(n: List, v: List) -> float:
    """get gini.

    Parameters
    ----------
    n : List
        list of numbers 1
    v : List
        size distribution

    Returns
    -------
    float
        gini index value
    """
    values = []
    for i in range(0, len(n)):
        for j in range(0, v[i]):
            values.append(n[i])
    n = len(values)
    assert n > 0, "Empty list of values"
    # Sort smallest to largest
    sortedValues = sorted(values)
    cumm = [0]
    for i in range(n):
        cumm.append(sum(sortedValues[0 : (i + 1)]))
    LorenzPoints = [[], []]
    # Some of all y values
    sumYs = 0
    # Robin Hood index max(x_i, y_i)
    robinHoodIdx = -1
    for i in range(1, n + 2):
        x = 100.0 * (i - 1) / n
        y = 100.0 * (cumm[i - 1] / float(cumm[n]))
        LorenzPoints[0].append(x)
        LorenzPoints[1].append(y)
        sumYs += y
        maxX_Y = x - y
        if maxX_Y > robinHoodIdx:
            robinHoodIdx = maxX_Y
    # Gini index
    giniIdx = 100 + (100 - 2 * sumYs) / n
    return giniIdx / 100


def get_network_statistics_per_chain(
    cluster_file: Path,
    sample_id: str,
    per_chain_repertoire_statistics_file: Path,
):
    """Summary

    Parameters
    ----------
    cluster_file : Path
        path to clustered file.
    sample_id : str
        name of sample.
    per_chain_repertoire_statistics_file : Path
        path to output statistics file.
    """
    create_file(per_chain_repertoire_statistics_file)
    fh = open(cluster_file, "r")
    cluster = Tree()
    index, sizesv, c_sizes = 0, [], {}
    total_v, total_reads = [], []
    sizesv, c_sizes = {}, {}
    chains_short = []
    t1 = 0
    n = 0
    for l in fh:
        if l[0] != "#":
            l = l.strip().split()
            id = l[2]
            chains, freq, id_short = (
                id.split("|")[1].split("_"),
                list(
                    map(
                        int,
                        id.split(READ_NUMBER_DIVISION)[1]
                        .split("|")[0]
                        .split("_"),
                    )
                ),
                id.split(READ_NUMBER_DIVISION)[0],
            )
            t1 = t1 + sum(freq)
            n = n + 1
            if len(chains_short) == 0:
                for i in range(0, len(chains)):
                    c = chains[i].split("*")[0]
                    if c not in chains_short:
                        chains_short.append(c)
            non_zero = [i for i in range(len(freq)) if freq[i] != 0]
            if len(total_v) == 0:
                total_v, total_reads = [0] * len(chains_short), [0] * len(
                    chains_short
                )
                for c in chains_short:
                    sizesv[c], c_sizes[c] = [], []
            for i in non_zero:
                c = chains[i].split("*")[0]
                cluster[c][l[1]][freq[i]][id_short].value = 1
                index = chains_short.index(c)
                total_v[index] = total_v[index] + 1
                sizesv[c] = sizesv[c] + [freq[i]]
                total_reads[index] = total_reads[index] + freq[i]
    fh.close()
    print(total_reads, t1, n)
    if t1 != sum(total_reads):
        print("ERROR IN COUNTING!!")
    out = (
        "#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\t"
        + "Cluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\n"
    )
    for c1 in chains_short:
        cluster_sizes_sub = []
        for clus in cluster[c1]:
            f = 0
            for f1 in cluster[c1][clus]:
                f = f + (f1 * len(cluster[c1][clus][f1]))
            cluster_sizes_sub = cluster_sizes_sub + [f]
        if len(cluster_sizes_sub) > 0:
            (vpoints, vvdf) = vdf(n=sizesv[c1])
            (cpoints, cvdf) = vdf(n=cluster_sizes_sub)
            vgini, cgini = gini_index(
                cpoints=cpoints,
                cvdf=cvdf,
                vpoints=vpoints,
                vvdf=vvdf,
            )
            max_pop, max_1_pop = cpoints[len(cpoints) - 1] * 100.0 / sum(
                sizesv[c1]
            ), cpoints[len(cpoints) - 2] * 100.0 / sum(sizesv[c1])
            out = (
                out
                + str(sample_id)
                + "\t"
                + c1
                + "\t"
                + str(sum(sizesv[c1]))
                + "\t"
                + str(len(sizesv[c1]))
                + "\t"
                + str(vgini)
                + "\t"
                + str(cgini)
                + "\t"
                + str(max_pop)
                + "\t"
                + str(max_1_pop)
                + "\n"
            )
    write_out(out, per_chain_repertoire_statistics_file)


def get_cluster_vertex_distributions(cluster_file: Path) -> Tuple:
    """Read clustered file and obtain stats

    Parameters
    ----------
    cluster_file : Path
        path to cluster file.

    Returns
    -------
    Tuple
        general statistics about cluster vertices,
    """
    fh = open(cluster_file, "r")
    cluster = Tree()
    (
        index,
        totalc,
        totalv,
        totalreads,
        sizesv,
        c_sizes,
        vertices_in_max_cluster,
    ) = (0, 0, 0, 0, [], {}, 0)
    for l in fh:
        index = index + 1
        if index > 1:
            l = l.strip()
            l = l.split()
            cluster[l[1]][l[2]].value = 1
            size = int(l[3])
            sizesv.append(size)
            totalv = totalv + 1
            totalreads = totalreads + size
            if int(l[1]) == 1:
                vertices_in_max_cluster = vertices_in_max_cluster + 1
            if l[1] in c_sizes:
                c_sizes[l[1]] = c_sizes[l[1]] + size
            else:
                c_sizes[l[1]] = size
    fh.close()
    sizes = []
    totalc = len(cluster)
    for c in cluster:
        sizes.append(len(cluster[c]))
    (cpoints, cvdf) = vdf(sizes)
    (vpoints, vvdf) = vdf(sizesv)
    return (
        cpoints,
        cvdf,
        vpoints,
        vvdf,
        totalc,
        totalv,
        totalreads,
        c_sizes,
        vertices_in_max_cluster,
    )


def proportional_measures(
    c_sizes: List, totalreads: int
) -> Tuple[float, float]:
    """Summary

    Parameters
    ----------
    c_sizes : List
        lsit of sizes
    totalreads : int
        total read count

    Returns
    -------
    Tuple[float, float]
        proportions
    """
    sizes = []
    for c in c_sizes:
        sizes.append((c, c_sizes[c]))
    s = sorted(sizes, key=itemgetter(1), reverse=True)
    (max_pop, max_1_pop) = (
        s[0][1] * 100.0 / totalreads,
        s[1][1] * 100.0 / totalreads,
    )
    return (max_pop, max_1_pop)


def print_distributions(points: List, cdf: List, file_out: Path):
    """Print distributions.

    Parameters
    ----------
    points : List
        list of sizes
    cdf : List
        cumulative distribution
    file_out : Path
        path to output
    """
    create_file(file_out)
    out = "#size\tfrequency of size\n"
    for i in range(0, len(points)):
        out = out + str(points[i]) + "\t" + str(cdf[i]) + "\n"
    write_out(out, file_out)


def get_network_statistics(
    cluster_file: Path,
    sample_id: str,
    network_statistics: Path,
    species: str,
    cluster_size_distribution: Path,
    vertex_size_distribution: Path,
):
    """Get network statistics

    Parameters
    ----------
    cluster_file : Path
        path to clustered file.
    sample_id : str
        name of sample.
    network_statistics : Path
        path to output statistics file.
    species : str
        organism type.
    cluster_size_distribution : Path
        path to output statistics file (cluster size distribution).
    vertex_size_distribution : Path
        path to output statistics file (vertex size distribution).
    """
    create_file(network_statistics)
    (
        cpoints,
        cvdf,
        vpoints,
        vvdf,
        totalc,
        totalv,
        totalreads,
        c_sizes,
        vertices_in_max_cluster,
    ) = get_cluster_vertex_distributions(cluster_file=cluster_file)
    print_distributions(
        points=vpoints, cdf=vvdf, file_out=vertex_size_distribution
    )
    print_distributions(
        points=cpoints, cdf=cvdf, file_out=cluster_size_distribution
    )
    vrenyi, crenyi = renyi_entropy(
        cpoints=cpoints,
        cvdf=cvdf,
        vpoints=vpoints,
        vvdf=vvdf,
        totalv=totalv,
        totalreads=totalreads,
    )
    vgini, cgini = gini_index(
        cpoints=cpoints,
        cvdf=cvdf,
        vpoints=vpoints,
        vvdf=vvdf,
    )
    max_pop, max_1_pop = proportional_measures(
        c_sizes=c_sizes, totalreads=totalreads
    )
    out = (
        "#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t"
        + "2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
    )
    out = (
        out
        + str(sample_id)
        + "\tOVERALL\t"
        + str(totalreads)
        + "\t"
        + str(totalv)
        + "\t"
        + str(vgini)
        + "\t"
        + str(cgini)
        + "\t"
        + str(max_pop)
        + "\t"
        + str(max_1_pop)
        + "\t"
        + str(vertices_in_max_cluster * 100.0 / totalv)
        + "\t"
        + str(vrenyi)
        + "\t"
        + str(crenyi)
        + "\t"
        + "IGH"  # always IGH in Clatworthy lab
        + "\t"
        + species
        + "\n"
    )
    write_out(out, network_statistics)


def main():
    """main function for step 4."""
    get_network_statistics(
        cluster_file=CLUST_ID_FILE,
        sample_id=SAMPLE_ID,
        network_statistics=NETSTATS,
        species=ORG,
        cluster_size_distribution=CLUSTER_SIZE_DIST,
        vertex_size_distribution=VERTEX_SIZE_DIST,
    )
    get_network_statistics_per_chain(
        cluster_file=CLUST_ID_FILE,
        sample_id=SAMPLE_ID,
        per_chain_repertoire_statistics_file=NETSTATS_PER_CHAIN,
    )


if __name__ == "__main__":
    main()
