#!/usr/bin/env python

# To do:
# - allow plot by genes


import os, sys

def read_matrix(filename):
    import arrayio

    return arrayio.read(filename)

def main():
    from optparse import OptionParser, OptionGroup
    import arrayio
    from genomicode import jmath
    from genomicode import pcafns
    from genomicode import parsefns
    from genomicode import colorfns

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
    
    # Parse the input arguments.
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("Please specify an infile and an outfile.")
    filename, outfile = args
    if not os.path.exists(filename):
        parser.error("I could not find file %s." % filename)

    num_genes = options.genes or 200
    K = 3  # number of dimensions

    MATRIX = read_matrix(filename)

    if options.log_transform:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)

    # Figure out the clusters for each of the samples.
    index2cluster = {}
    for clust_i, s in enumerate(options.cluster):
        ranges = parsefns.parse_ranges(s)
        for s, e in ranges:
            for i in range(s-1, e):
                assert i < MATRIX.ncol(), "Index %d out of range" % i
                assert i not in index2cluster, \
                       "Index %d in multiple clusters" % i
                index2cluster[i] = clust_i
    cluster = [len(options.cluster)] * MATRIX.ncol()
    for i, g in index2cluster.iteritems():
        cluster[i] = g

    # Select a subset of the genes.
    I = pcafns.select_genes_var(MATRIX._X, num_genes)
    MATRIX = MATRIX.matrix(I, None)

    # Calculate the principal components and plot them.
    principal_components = pcafns.svd_project_cols(MATRIX._X, K)
    X = [x[0] for x in principal_components]
    Y = [x[1] for x in principal_components]
    color = pcafns.choose_colors(cluster)
    pcafns.plot_scatter(X, Y, outfile, group=cluster, color=color)

    # Write out the principal components.
    assert len(principal_components) == len(cluster)
    x = ["PC%d" % i for i in range(K)]
    x = ["Index", "Sample", "Cluster", "Color"] + x
    print "\t".join(x)
    for i in range(len(principal_components)):
        x = MATRIX.col_names(arrayio.COL_ID)[i]
        c = ""
        if color:
            c = colorfns.rgb2hex(color[i])
        x = [i+1, x, cluster[i], c] + principal_components[i]
        print "\t".join(map(str, x))


if __name__ == '__main__':
    main()

