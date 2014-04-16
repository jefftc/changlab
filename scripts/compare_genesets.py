#!/usr/bin/env python

import os
import sys

def _parse_file_gs(geneset):
    # Parse a geneset specified by the user.  geneset is in the format
    # of <filename>[,<geneset>,<geneset>,...].  Return a tuple of
    # <filename>, list of <geneset> (or empty list).
    # XXX what happens if this is an empty list?
    x = geneset.split(",")
    assert len(x) >= 1
    filename, genesets = x[0], x[1:]
    return filename, genesets


def main():
    import argparse

    from genomicode import genesetlib

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--operation", choices=["merge", "intersect"])
    parser.add_argument(
        "--geneset", default=[], action="append",
        help="Format: <gmx/gmt_file>,<geneset>")
    parser.add_argument("--geneset_name")
    parser.add_argument("--format", default="gmt", choices=["gmt", "gmx"])
    
    args = parser.parse_args()
    assert len(args.geneset) >= 2, "At least two gene sets."

    # Read the gene sets.
    name2genes = {}
    for x in args.geneset:
        x = _parse_file_gs(x)
        filename, genesets = x
        assert os.path.exists(filename), "File not found: %s" % filename
        assert len(genesets) == 1
        geneset = genesets[0]

        genes = genesetlib.read_genes(filename, geneset)
        name2genes[geneset] = genes

    # Come up with a geneset name.
    geneset_name = args.geneset_name
    if not geneset_name:
        geneset_name = args.operation
    
    genesets = []
    if args.operation == "merge":
        all_genes = []
        for name, genes in name2genes.iteritems():
            all_genes.extend(genes)
        all_genes = sorted({}.fromkeys(all_genes))
        x = genesetlib.GeneSet(geneset_name, "na", all_genes)
        genesets.append(x)
    elif args.operation == "intersect":
        gene2count = {}
        for name, genes in name2genes.iteritems():
            genes = {}.fromkeys(genes)
            for x in genes:
                gene2count[x] = gene2count.get(x, 0) + 1

        intersect_genes = {}
        for gene, count in gene2count.iteritems():
            if count == len(name2genes):
                intersect_genes[gene] = 1
        intersect_genes = sorted(intersect_genes)

        name2intersect = {}
        for name, genes in name2genes.iteritems():
            genes = {}.fromkeys(genes)
            genes = [x for x in genes if x not in intersect_genes]
            name2intersect[name] = genes

        x = genesetlib.GeneSet(geneset_name, "na", sorted(intersect_genes))
        genesets.append(x)
        for name, genes in name2intersect.iteritems():
            if not genes:
                continue
            x = genesetlib.GeneSet("%s_LEFTOVER" % name, "na", sorted(genes))
            genesets.append(x)
    else:
        raise NotImplementedError, args.operation

    if args.format == "gmt":
        genesetlib.write_gmt(sys.stdout, genesets)
    elif args.format == "gmx":
        genesetlib.write_gmx(sys.stdout, genesets)
    else:
        raise NotImplementedError, args.format


if __name__ == '__main__':
    main()
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals()

