#!/usr/bin/env python

# To do:
# - allow plot by genes


import os, sys

def read_matrix(filename):
    import arrayio

    return arrayio.read(filename)

def _parse_cluster(options_cluster, MATRIX):
    # Return a vector of clusters, where each cluster is an integer
    # from 0 to K-1.  K is the total number of clusters.  The length
    # of the vector should be the same as the number of samples in the
    # matrix.
    from genomicode import parselib

    index2cluster = {}
    for clust_i, s in enumerate(options_cluster):
        ranges = parselib.parse_ranges(s)
        for s, e in ranges:
            for i in range(s-1, e):
                assert i < MATRIX.ncol(), "Index %d out of range" % i
                assert i not in index2cluster, \
                       "Index %d in multiple clusters" % i
                index2cluster[i] = clust_i
    #cluster = [len(options_cluster)] * MATRIX.ncol()
    cluster = [None] * MATRIX.ncol()
    for i, g in index2cluster.iteritems():
        cluster[i] = g
    return cluster

def _parse_cluster_file(cluster_file, MATRIX):
    # Return a vector of clusters, where each cluster is an integer
    # from 0 to K-1.  K is the total number of clusters.  The length
    # of the vector should be the same as the number of samples in the
    # matrix.
    from genomicode import clusterio
    
    id2cluster = {}
    for (id, cluster) in clusterio.read_kgg_file(cluster_file):
        id2cluster[id] = cluster

    # Figure out which row header matches the IDs.
    header2numids = {}  # header -> number of IDs matched.
    header = num_ids = None
    for cn in MATRIX.col_names():
        x = MATRIX.col_names(cn)
        col_ids = {}.fromkeys(x)
        x = [x for x in id2cluster if x in col_ids]
        if header is None or len(x) > num_ids:
            header, num_ids = cn, len(x)

    index2cluster = {}
    if header is not None:
        for i, name in enumerate(MATRIX.col_names(header)):
            cluster = id2cluster.get(name)
            index2cluster[i] = cluster
            
    cluster = [None] * MATRIX.ncol()
    for i, g in index2cluster.iteritems():
        cluster[i] = g

    # Make sure the clusters are 0-based.
    clean = [x for x in cluster if x is not None]
    if clean and min(clean) > 0:
        for i in range(len(cluster)):
            if cluster[i] is None:
                continue
            cluster[i] = cluster[i] - min(clean)
    return cluster

def main():
    from optparse import OptionParser, OptionGroup
    import arrayio
    from genomicode import jmath
    from genomicode import pcalib
    from genomicode import colorlib

    usage = "usage: %prog [options] filename outfile"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "-l", "--log_transform", dest="log_transform", default=False,
        action="store_true",
        help="Log transform the data first.")
    parser.add_option(
        "-g", "--genes", dest="genes", default=None, type="int",
        help="Number of genes to use.")
    parser.add_option(
        "-c", "--cluster", dest="cluster", default=[], action="append",
        help="Group samples into a cluster (e.g. -c 1-5); 1-based.")
    parser.add_option(
        "", "--cluster_file", dest="cluster_file", default=None, 
        help="A KGG format file of the clusters for the samples.  "
        "Clusters in this file can be 0-based or 1-based.")
    parser.add_option(
        "", "--label", dest="label", default=False, action="store_true",
        help="Label the samples.")
    parser.add_option(
        "", "--title", dest="title", default=None, type="str",
        help="Put a title on the plot.")
    parser.add_option(
        "-v", "--verbose", dest="verbose", default=False, action="store_true",
        help="")
    
    # Parse the input arguments.
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Please specify an infile and an outfile.")
    elif len(args) > 2:
        parser.error("Too many input parameters (%d)." % len(args))
    filename, outfile = args
    if not os.path.exists(filename):
        parser.error("I could not find file %s." % filename)

    num_genes = options.genes
    #K = 10  # number of dimensions

    MATRIX = read_matrix(filename)
    if options.log_transform:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)

    cluster = None
    if options.cluster and options.cluster_file:
        parser.error("Cannot specify clusters and a cluster file.")
    if options.cluster:
        cluster = _parse_cluster(options.cluster, MATRIX)
    if options.cluster_file:
        if not os.path.exists(options.cluster_file):
            parser.error(
                "I could not find cluster file: %s" % options.cluster_file)
        cluster = _parse_cluster_file(options.cluster_file, MATRIX)

    # Select a subset of the genes.
    if num_genes:
        I = pcalib.select_genes_var(MATRIX._X, num_genes)
        MATRIX = MATRIX.matrix(I, None)

    # Calculate the principal components and plot them.
    K = min(MATRIX.nrow(), MATRIX.ncol())
    principal_components, perc_var = pcalib.svd_project_cols(MATRIX._X, K)
    X = [x[0] for x in principal_components]
    Y = [x[1] for x in principal_components]
    color = None
    if cluster is not None:
        color = pcalib.choose_colors(cluster)
    LABEL = None
    if options.label:
        LABEL = MATRIX.col_names(arrayio.COL_ID)
    assert not LABEL or len(LABEL) == len(X), "%d %d" % (len(X), len(LABEL))
    pcalib.plot_scatter(
        X, Y, outfile, group=cluster, color=color, title=options.title,
        label=LABEL)

    if options.verbose:
        # Write out the principal components.
        assert cluster is None or len(cluster) == len(principal_components)
        x = ["PC%02d (%.2f%%)" % (i, 100*perc_var[i]) for i in range(K)]
        header = ["Index", "Sample", "Cluster", "Color"] + x
        print "\t".join(header)
        for i in range(len(principal_components)):
            x = MATRIX.col_names(arrayio.COL_ID)[i]
            c = ""
            if color:
                c = colorlib.rgb2hex(color[i])
            clust = ""
            if cluster is not None:
                clust = cluster[i]
            x = [i+1, x, clust, c] + principal_components[i]
            assert len(x) == len(header)
            print "\t".join(map(str, x))


if __name__ == '__main__':
    main()

