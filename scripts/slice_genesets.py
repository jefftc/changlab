#!/usr/bin/env python

import os
import sys

def merge_genesets(genesets, geneset_name):
    from genomicode import genesetlib
    if not geneset_name:
        return genesets
    genesets = genesets[:]

    all_genes = []
    for gs in genesets:
        all_genes.extend(gs.genes)
    all_genes = sorted({}.fromkeys(all_genes))
    x = genesetlib.GeneSet(geneset_name, "na", all_genes)
    genesets.append(x)
    return genesets


def intersect_genesets(genesets, geneset_name):
    from genomicode import genesetlib
    if not geneset_name:
        return genesets
    genesets = genesets[:]

    gene2count = {}
    for gs in genesets:
        genes = {}.fromkeys(gs.genes)
        for x in genes:
            gene2count[x] = gene2count.get(x, 0) + 1

    intersect_genes = {}
    for gene, count in gene2count.iteritems():
        if count == len(genesets):
            intersect_genes[gene] = 1
    intersect_genes = sorted(intersect_genes)

    name2intersect = {}
    for gs in genesets:
        genes = {}.fromkeys(gs.genes)
        genes = [x for x in genes if x not in intersect_genes]
        name2intersect[name] = genes

    x = genesetlib.GeneSet(geneset_name, "na", sorted(intersect_genes))
    genesets.append(x)
    for name, genes in name2intersect.iteritems():
        if not genes:
            continue
        x = genesetlib.GeneSet("%s_LEFTOVER" % name, "na", sorted(genes))
        genesets.append(x)
    
    return genesets


def clean_excel(genesets, clean):
    from genomicode import genesetlib
    if not clean:
        return genesets
    genesets = genesets[:]

    TO_FIX = {
        "3/7/07 0:00" : "LDOC1",  # (was MAR7?)
        "9/9/14 0:00" : "SEPT9",
        }
    for gs in genesets:
        gs.genes = [TO_FIX.get(x, x) for x in gs.genes]
    return genesets
    

def append_to_name(genesets, arg):
    if not arg:
        return genesets
    from genomicode import genesetlib
    genesets = genesets[:]

    for gs in genesets:
        gs.name = gs.name + arg
    return genesets


def rename(genesets, args):
    # args is a list of strings in format of: <from>,<to>.
    if not args:
        return genesets
    from genomicode import genesetlib

    rename_all = []  # list of (from_str, to_str)
    for rename_str in args:
        x = rename_str.split(",")
        if len(x) > 2 or len(x) == 1:
            x = rename_str.split(";")
        assert len(x) == 2, "format is not <from>,<to>: %s" % rename_str
        from_str, to_str = x
        rename_all.append((from_str, to_str))

    genesets = genesets[:]
    for gs in genesets:
        for from_str, to_str in rename_all:
            if gs.name == from_str:
                gs.name = to_str
    return genesets


def nodup(genesets, arg):
    if not arg:
        return genesets
    genesets = genesets[:]

    for gs in genesets:
        seen = {}
        genes = []
        for g in gs.genes:
            if g in seen:
                continue
            seen[g] = 1
            genes.append(g)
        gs.genes = genes
    return genesets


def sort_genes(genesets, arg):
    from genomicode import genesetlib
    from genomicode import sortlib
    if not arg:
        return genesets
    genesets = genesets[:]

    for gs in genesets:
        gs.genes = sortlib.sort_natural(gs.genes)
    return genesets


def sort_genesets(genesets, arg):
    from genomicode import sortlib
    if not arg:
        return genesets

    x = [x.name for x in genesets]
    O = sortlib.order_natural(x)
    genesets = [genesets[i] for i in O]
    return genesets


def _parse_file_gs(geneset):
    # Parse a geneset specified by the user.  geneset is in the format
    # of <filename>[,<geneset>,<geneset>,...].  Return a tuple of
    # <filename>, list of <geneset> (or empty list).
    x = geneset.split(",")
    assert len(x) >= 1
    filename, genesets = x[0], x[1:]
    return filename, genesets


def main():
    import argparse
    from genomicode import genesetlib

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "geneset", nargs="+", 
        help="Format: <gmx/gmt_file>[,<geneset>,<geneset>].  "
        "If <geneset> is not given, then I will use every gene set "
        "in this file.")
    DEFAULT_OUTPUT_FORMAT = "gmt"
    parser.add_argument(
        "--format", default=DEFAULT_OUTPUT_FORMAT, choices=["gmt", "gmx"],
        help="Output format.  Default %s." % DEFAULT_OUTPUT_FORMAT)

    group = parser.add_argument_group(title='Operations')
    group.add_argument(
        "--merge",
        help="Merge gene sets.  The argument should be the name of the "
        "new merged gene set.  Format: <geneset_name>.")
    group.add_argument(
        "--intersect",
        help="Merge gene sets.  The argument should be the name of the "
        "intersection.  Format: <geneset_name>.")
    group.add_argument(
        "--clean_excel", action="store_true",
        help="Clean up gene names that Excel messes up (like OCT1).")
    group.add_argument(
        "--append_to_name", 
        help="Append to the name of these gene sets.  "
        "Format: <text_to_append>.")
    group.add_argument(
        "--rename",  action="append",
        help="Change the name of a gene set.  Format: <from>,<to>.  "
        "<from> will be replaced with <to>.  "
        "If there are already commas in the header names, can use ; instead.  "
        "(MULTI)")
    group.add_argument(
        "--nodup", action="store_true", help="Remove duplicate genes.")
    group.add_argument("--sort_genes", action="store_true")
    group.add_argument("--sort_genesets", action="store_true")

    args = parser.parse_args()
    assert args.geneset
    #assert not (args.geneset and args.all_genesets)
    #assert len(args.geneset) >= 2, "At least two gene sets."

    # Read the gene sets.
    genesets = []
    all_names = []
    for x in args.geneset:
        x = _parse_file_gs(x)
        filename, names = x
        assert os.path.exists(filename), "File not found: %s" % filename

        all_names.extend(names)
        
        for x in genesetlib.read_genesets(filename):
            name, desc, genes = x
            if not names or name in names:
                x = genesetlib.GeneSet(name, desc, genes)
                genesets.append(x)
    assert genesets, "No gene sets"

    all_gs = [x.name for x in genesets]
    x = all_names
    x = [x for x in x if x not in all_gs]
    missing = x
    assert not missing, "Missing gene sets: %s" % missing

    genesets = merge_genesets(genesets, args.merge)
    genesets = intersect_genesets(genesets, args.intersect)
    genesets = clean_excel(genesets, args.clean_excel)
    genesets = append_to_name(genesets, args.append_to_name)
    genesets = rename(genesets, args.rename)
    genesets = nodup(genesets, args.nodup)
    genesets = sort_genes(genesets, args.sort_genes)
    genesets = sort_genesets(genesets, args.sort_genesets)


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
