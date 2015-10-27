#!/usr/bin/env python

# To do:
# - allow plot by genes


import os, sys

def read_matrix(filename, num_header_cols=None):
    import arrayio

    return arrayio.read(filename, hcols=num_header_cols)

def _parse_cluster(options_cluster, indexes_include_headers, MATRIX):
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
                if indexes_include_headers:
                    i -= len(MATRIX._row_names)
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
    from genomicode import prismlib

    usage = "usage: %prog [options] filename outfile.png"
    parser = OptionParser(usage=usage, version="%prog 01")

    #parser.add_option(
    #    "-l", "--log_transform", default=False,
    #    action="store_true",
    #    help="Log transform the data first.")
    
    parser.add_option(
        "--num_header_cols", type=int,
        help="This number of columns are headers.  If not given, will guess.")
    parser.add_option(
        "-g", "--genes", default=None, type="int",
        help="Number of genes to use.")
    parser.add_option(
        "--prism_file", 
        help="Write the results out to a prism-formatted file.")
    parser.add_option(
        "-v", "--verbose", default=False, action="store_true",
        help="")

    group = OptionGroup(parser, "Clustering")
    parser.add_option_group(group)
    group.add_option(
        "-c", "--cluster", default=[], action="append",
        help="Group samples into a cluster (e.g. -c 1-5); 1-based.")
    group.add_option(
        "--indexes_include_headers", "--iih", action="store_true",
        help="If not given (default), then index 1 is the first column "
        "with data.  If given, then index 1 is the very first column "
        "in the file, including the headers.")
    group.add_option(
        "--cluster_file", 
        help="A KGG format file of the clusters for the samples.  "
        "Clusters in this file can be 0-based or 1-based.")
    

    group = OptionGroup(parser, "Visualization")
    parser.add_option_group(group)
    group.add_option(
        "--title", help="Put a title on the plot.")
    group.add_option(
        "--width", default=None, type="int",
        help="Width (in pixels) of the plot.")
    group.add_option(
        "--label", default=False, action="store_true",
        help="Label the samples.")
    group.add_option(
        "--scale_label", type=float, default=1.0, 
        help="Scale the size of the labels.")
    
    # Parse the input arguments.
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Please specify an infile and an outfile.")
    elif len(args) > 2:
        parser.error("Too many input parameters (%d)." % len(args))
    filename, outfile = args
    if not os.path.exists(filename):
        parser.error("I could not find file %s." % filename)
    if options.num_header_cols is not None:
        assert options.num_header_cols > 0 and options.num_header_cols < 100
    if options.width is not None:
        assert options.width > 10, "too small"
        assert options.width < 4096*16, "width too big"
    assert options.scale_label > 0.01 and options.scale_label < 100
    options.log_transform = False

    num_genes = options.genes
    #K = 10  # number of dimensions

    MATRIX = read_matrix(filename, options.num_header_cols)
    if options.log_transform:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    cluster = None
    if options.cluster and options.cluster_file:
        parser.error("Cannot specify clusters and a cluster file.")
    if options.cluster:
        cluster = _parse_cluster(
            options.cluster, options.indexes_include_headers, MATRIX)
    if options.cluster_file:
        if not os.path.exists(options.cluster_file):
            parser.error(
                "I could not find cluster file: %s" % options.cluster_file)
        cluster = _parse_cluster_file(options.cluster_file, MATRIX)

    # Select a subset of the genes.
    if num_genes:
        assert MATRIX.ncol() > 1, "Not enough samples to select genes."
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
    height = width = None
    if options.width is not None:
        height, width = int(options.width*0.75), options.width
    pcalib.plot_scatter(
        X, Y, outfile, group=cluster, color=color, title=options.title,
        label=LABEL, xlabel=False, ylabel=False,
        scale_label=options.scale_label, height=height, width=width)

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

    if options.prism_file:
        # Write out as prism format.
        num_series = 1
        if cluster:
            num_series = max(cluster) + 1
        names = ["CLUSTER-%d" % (i+1) for i in range(num_series)]
        DATA = {}
        rownames = {}
        for i in range(num_series):
            xy = []
            n = []
            for j in range(len(principal_components)):
                if cluster and cluster[j] != i:
                    continue
                x = principal_components[j][0]
                y = principal_components[j][1]
                xy.append([x, y])
                n.append(MATRIX.col_names(arrayio.COL_ID)[j])
            if xy:
                DATA[names[i]] = xy
                rownames[names[i]] = n

        prismlib.write_scatterplot(options.prism_file, DATA, rownames)


if __name__ == '__main__':
    main()

