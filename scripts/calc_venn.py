#!/usr/bin/env python

def main():
    import os
    import argparse
    import itertools
    from genomicode import genesetlib

    # Option for case insensitive comparison?
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "geneset_file", help="File in GMX or GMT format.")
    parser.add_argument(
        "-o", dest="outfile", help="Save the intersection to this "
        "(GMX or GMT) file.")

    args = parser.parse_args()
    assert os.path.exists(args.geneset_file), \
           "File not found: %s" % args.geneset_file

    name2genes = {}
    all_names = []  # preserve order of gene sets
    for x in genesetlib.read_genesets(args.geneset_file):
        name, desc, genes = x
        name2genes[name] = genes
        all_names.append(name)

    # Count all pairwise intersections.
    pair2common = {}  # (name1, name2) -> common genes
    for x in itertools.product(all_names, all_names):
        name1, name2 = x
        genes1 = name2genes[name1]
        genes2 = name2genes[name2]
        common = set(genes1).intersection(genes2)
        pair2common[(name1, name2)] = sorted(common)

    # Write out the matrix.
    header = ["Count"] + all_names
    print "\t".join(header)
    for name1 in all_names:
        x = [len(pair2common[(name1, name2)]) for name2 in all_names]
        x = [name1] + x
        assert len(x) == len(header)
        print "\t".join(map(str, x))

    # Write to an output file.
    if not args.outfile:
        return
    genesets = []  # list of GeneSet objects
    for i in range(len(all_names)-1):
        for j in range(i, len(all_names)):
            name1 = all_names[i]
            name2 = all_names[j]
            genes = pair2common[(name1, name2)]
            name = "%s_%s" % (name1, name2)
            x = genesetlib.GeneSet(name, "na", genes)
            genesets.append(x)

    write_fn = genesetlib.write_gmx
    if args.outfile.lower().endswith(".gmt"):
        write_fn = genesetlib.write_gmt
    write_fn(args.outfile, genesets)
            

if __name__ == '__main__':
    main()
