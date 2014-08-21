#!/usr/bin/env python

def main():
    import os
    import argparse
    import itertools
    from genomicode import genesetlib

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "geneset_file", help="File in GMX or GMT format.")
    # Case sensitive?

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
    pair2count = {}  # (name1, name2) -> number of common genes
    for x in itertools.product(all_names, all_names):
        name1, name2 = x
        genes1 = name2genes[name1]
        genes2 = name2genes[name2]
        common = set(genes1).intersection(genes2)
        pair2count[(name1, name2)] = len(common)

    # Write out the matrix.
    header = ["Count"] + all_names
    print "\t".join(header)
    for name1 in all_names:
        x = [pair2count[(name1, name2)] for name2 in all_names]
        x = [name1] + x
        assert len(x) == len(header)
        print "\t".join(map(str, x))
            

if __name__ == '__main__':
    main()
