#!/usr/bin/env python

# Functions:

import os
import sys


def read_matrix(filename):
    import arrayio

    # Read the files.
    assert os.path.exists(filename), \
        "I could not find the file: %s" % filename
    fmt_module = arrayio.choose_format(filename)
    assert fmt_module, \
        "I could not figure out the format of file: %s" % filename
    x = fmt_module.read(filename)
    return x


def write_matrix(filename, matrix):
    import arrayio
    arrayio.write(matrix, open(filename, 'w'))


def write_annot(filename, name2annots):
    from genomicode import jmath
    matrix = []
    for name in sorted(name2annots):
        annots = name2annots[name]
        x = [name] + annots
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)

    handle = open(filename, 'w')
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


def read_annot(filename):
    from genomicode import genesetlib

    name2annots = {}
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, annots = x
        name2annots[name] = annots
    return name2annots


def get_sample_names(matrix):
    import arrayio
    
    name = arrayio.COL_ID
    if name not in matrix._col_names:
        name = matrix._synonyms[name]
    assert name in matrix._col_names, "I can not find the sample names."
    x = matrix.col_names(name)
    return x


def find_best_annots(name2annots, match_annots):
    # Return list of (name, matched_annotations)
    match_annots_s = set(match_annots)
    all_matches = []  # list of (num_matches, name, matches)
    for name, annots in name2annots.iteritems():
        matches = list(match_annots_s.intersection(annots))
        x = len(matches), name, matches
        all_matches.append(x)
    all_matches = sorted(all_matches)
    x = all_matches[-1]
    x, name, matches = x
    return name, matches
    

def find_common_samples(matrix_data):
    # First, get the common samples from all the expression data.
    num_express = 0
    sample2count = {}
    for x in matrix_data:
        infile, outfile, is_express_file, matrix = x
        if not is_express_file:
            continue
        num_express += 1
        for x in get_sample_names(matrix):
            sample2count[x] = sample2count.get(x, 0) + 1
    if not num_express:
        raise NotImplementedError
            
    common = [s for (s, c) in sample2count.iteritems() if c == num_express]
    if not common:
        return []

    # Now, match up with the annotation files.
    num_express = 0
    sample2count = {}
    for x in matrix_data:
        infile, outfile, is_express_file, name2annots = x
        if is_express_file:
            continue
        name, annots = find_best_annots(name2annots, common)
        assert annots, "%s has no matching annotations." % infile
        common = annots

    return sorted(common)


def align_express(matrix, samples):
    names = get_sample_names(matrix)
    I = []
    for n in samples:
        assert n in names
        i = names.index(n)
        I.append(i)
    x = matrix.matrix(None, I)
    return x


def align_annot(name2annots, samples):
    name, x = find_best_annots(name2annots, samples)
    annots = name2annots[name]
    I = []
    for n in samples:
        assert n in annots
        i = annots.index(n)
        I.append(i)
    name2annots_new = {}
    for name, annots in name2annots.iteritems():
        x = [annots[i] for i in I]
        name2annots_new[name] = x
    return name2annots_new


def align_matrices(matrix_data, samples):
    new_matrix_data = []
    for x in matrix_data:
        infile, outfile, is_express_file, data = x
        if is_express_file:
            data_new = align_express(data, samples)
        else:
            data_new = align_annot(data, samples)
        x = infile, outfile, is_express_file, data_new
        new_matrix_data.append(x)
    return new_matrix_data


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Align a set of matrices.")
    parser.add_argument("outfile", nargs="+")
    
    parser.add_argument(
        "--clobber", default=False, action="store_true",
        help="Overwrite output files, if they already exist.")

    parser.add_argument(
        "--express_file", default=[], action="append", help="")
    parser.add_argument(
        "--annot_file", default=[], action="append", help="")
    
    args = parser.parse_args()
    assert len(args.outfile) == len(args.express_file) + len(args.annot_file)
    for x in args.express_file + args.annot_file:
        assert os.path.exists(x), "I could not find file: %s" % x
    for x in args.outfile:
        assert args.clobber or not os.path.exists(x), "File exists: %s" % x

    # Align the outfiles to the expression and annotation files.
    express_file = args.express_file[:]
    annot_file = args.annot_file[:]
    outfile = args.outfile[:]
    matrix_data = []  # list of (infile, outfile, is_express_file)
    for arg in sys.argv:
        if arg not in ["--express_file", "--annot_file"]:
            continue
        assert outfile
        if arg == "--express_file":
            assert express_file
            x = express_file.pop(0), outfile.pop(0), True
        else:
            assert annot_file
            x = annot_file.pop(0), outfile.pop(0), False
        matrix_data.append(x)
    assert not express_file
    assert not annot_file
    assert not outfile

    # Read each of the files.
    new_matrix_data = []  # list of (infile, outfile, is_express_file, data)
    for x in matrix_data:
        infile, outfile, is_express_file = x
        if is_express_file:
            data = read_matrix(infile)
        else:
            data = read_annot(infile)
        x = infile, outfile, is_express_file, data
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    # Find the intersection of each file.
    common_samples = find_common_samples(matrix_data)
    assert common_samples, "No common samples found."

    # Align each of the matrices.
    matrix_data = align_matrices(matrix_data, common_samples)

    # Write out each of the matrices.
    for x in matrix_data:
        infile, outfile, is_express_file, data = x
        if is_express_file:
            write_matrix(outfile, data)
        else:
            write_annot(outfile, data)
    

if __name__ == '__main__':
    main()
