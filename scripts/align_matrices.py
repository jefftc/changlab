#!/usr/bin/env python

# Functions:
# list_all_samples
# list_common_samples
# peek_samples_hint
# 
# align_matrices        Align a set of matrices.
# align_matrix          Align either expression or annotation matrices.
# align_express
# align_annot
#
# get_samples           Get the samples from an express or annot matrix.
# get_express_samples
# get_annot_samples
# 
# cmp_sample
# find_sample
# intersect_samples
#
# write_matrix
# read_express
# write_express
# read_annot
# write_annot

import os
import sys


class AnnotationMatrix:
    def __init__(self, name2annots, name_order):
        self.name2annots = name2annots.copy()
        self.name_order = name_order[:]


def list_all_samples(matrix_data, first_annot_header, case_insensitive):
    assert matrix_data

    samples_hint = peek_samples_hint(
        matrix_data, first_annot_header, case_insensitive)

    # Get the samples that occur in any of the files.  Preserve the
    # order of the samples.
    all_samples = []
    for x in matrix_data:
        infile, outfile, matrix = x
        samples = get_samples(matrix, samples_hint, case_insensitive)
        all_samples.extend(samples)

    # Get rid of duplicate samples.
    duplicated = [False] * len(all_samples)
    all_samples_cmp = all_samples
    if case_insensitive:
        all_samples_cmp = [x.upper() for x in all_samples]
    seen = {}
    for i, s in enumerate(all_samples_cmp):
        if s in seen:
            duplicated[i] = True
        seen[s] = 1

    all_samples = [
        all_samples[i] for i in range(len(all_samples)) if not duplicated[i]]
    return all_samples


def list_common_samples(matrix_data, first_annot_header, case_insensitive):
    assert matrix_data

    # Initialize the common samples.
    common = peek_samples_hint(
        matrix_data, first_annot_header, case_insensitive)
    assert common

    # Get the samples that are in all of the data sets.
    for x in matrix_data:
        infile, outfile, matrix = x
        samples = get_samples(matrix, common, case_insensitive)
        common = intersect_samples(common, samples, case_insensitive)

    # Finally, order the annotations according to the first file given.
    assert matrix_data
    infile, outfile, matrix = matrix_data[0]
    samples = get_samples(matrix, common, case_insensitive)
    samples = intersect_samples(samples, common, case_insensitive)
    
    return samples


def peek_samples_hint(matrix_data, first_annot_header, case_insensitive):
    # Figure out what the samples look like.  Look for an expression
    # file first.  If an expression file is not found, use the
    # annotation files.
    samples_hint = None
    x = [x for x in matrix_data if not isinstance(x[2], AnnotationMatrix)]
    if x:
        # Found an expression file.
        infile, outfile, matrix = x[0]
        samples_hint = get_express_samples(matrix)
    else:
        # No expression files.
        assert first_annot_header, "--first_annot_header not provided."
        infile, outfile, matrix = matrix_data[0]
        assert first_annot_header in matrix.name2annots
        samples_hint = matrix.name2annots[first_annot_header]
    assert samples_hint
    return samples_hint


def align_matrices(matrix_data, samples, case_insensitive, null_string):
    new_matrix_data = []
    for x in matrix_data:
        infile, outfile, matrix = x
        aligned_matrix = align_matrix(
            matrix, samples, case_insensitive, null_string)
        x = infile, outfile, aligned_matrix
        new_matrix_data.append(x)
    return new_matrix_data


def align_matrix(matrix, samples, case_insensitive, null_string):
    if isinstance(matrix, AnnotationMatrix):
        aligned_matrix = align_annot(
            matrix, samples, case_insensitive, null_string)
    else:
        aligned_matrix = align_express(
            matrix, samples, case_insensitive, null_string)
    return aligned_matrix


def align_express(matrix, samples, case_insensitive, null_string):
    import arrayio
    
    names = get_express_samples(matrix)
    I = []
    for n in samples:
        i = find_sample(names, n, case_insensitive)
        if i == -1:
            i = None
        I.append(i)

    # The Matrix class requires actually indexes.  Give it index 0,
    # and then remove those lines myself later.
    I_index = []
    for i in I:
        if i != None:
            I_index.append(i)
        else:
            I_index.append(0)
    assert len(I_index) == len(I)
    matrix_aligned = matrix.matrix(None, I_index)

    # Fix the None indexes in the matrix.
    X = matrix_aligned._X
    assert len(I) == len(X[0])
    for i_new, i_old in enumerate(I):
        if i_old != None:
            continue
        # Fill column i with blanks.
        for j in range(len(X)):
            X[j][i_new] = ""

    # Fix the annotations.
    header = arrayio.COL_ID
    if header not in matrix._col_names:
        header = matrix._synonyms[header]
    assert header in matrix._col_names, "I can not find the sample names."
    
    # Fix the sample names.
    for name in matrix_aligned._col_names:
        annots = matrix_aligned._col_names[name]
        annots_new = []
        for i_new, i_old in enumerate(I):
            if i_old != None:
                annots_new.append(annots[i_new])
            elif name == header:
                annots_new.append(samples[i_new])
            else:
                annots_new.append(null_string)
        matrix_aligned._col_names[name] = annots_new
    
    return matrix_aligned


def align_annot(matrix, samples, case_insensitive, null_string):
    names = get_annot_samples(matrix, samples, case_insensitive)

    # Figure out which is the header for these samples.
    header = None
    for h, n in matrix.name2annots.iteritems():
        if n == names:
            header = h
    assert header is not None
    
    I = []
    for n in samples:
        i = find_sample(names, n, case_insensitive)
        if i == -1:
            i = None
        I.append(i)
    assert len(I) == len(samples)
        
    name2annots_new = {}
    for name, annots in matrix.name2annots.iteritems():
        annots_new = []
        for i, i_annot in enumerate(I):
            if i_annot != None:
                annots_new.append(annots[i_annot])
            elif name == header:
                annots_new.append(samples[i])
            else:
                annots_new.append(null_string)
        name2annots_new[name] = annots_new
    return AnnotationMatrix(name2annots_new, matrix.name_order)


def get_samples(matrix, samples_hint, case_insensitive):
    # Get the samples from the matrix.  Since in principle anything in
    # the annotation file can be a sample, need to give it a hint of
    # what the samples look like.
    if isinstance(matrix, AnnotationMatrix):
        return get_annot_samples(matrix, samples_hint, case_insensitive)
    return get_express_samples(matrix)


def get_express_samples(matrix):
    import arrayio
    
    name = arrayio.COL_ID
    if name not in matrix._col_names:
        name = matrix._synonyms[name]
    assert name in matrix._col_names, "I can not find the sample names."
    x = matrix.col_names(name)
    return x


def get_annot_samples(matrix, samples_hint, case_insensitive):
    all_matches = []  # list of (num_matches, name, matches)
    for name, annots in matrix.name2annots.iteritems():
        matches = intersect_samples(annots, samples_hint, case_insensitive)
        x = len(matches), name, matches
        all_matches.append(x)
    assert all_matches
    all_matches = sorted(all_matches)
    x = all_matches[-1]
    x, name, x = x
    return matrix.name2annots[name]


def cmp_sample(x, y, case_insensitive):
    if case_insensitive:
        x, y = x.upper(), y.upper()
    return x == y


## def index_sample(sample_list, sample, case_insensitive):
##     i = find_sample(sample_list, sample, case_insensitive)
##     assert i >= 0, "Missing: %s" % sample


def find_sample(sample_list, sample, case_insensitive):
    # Return the index of this sample or -1.
    sample_list_cmp = sample_list
    sample_cmp = sample
    if case_insensitive:
        sample_list_cmp = [x.upper() for x in sample_list]
        sample_cmp = sample.upper()
    #I = []
    for i, x in enumerate(sample_list_cmp):
        if x == sample_cmp:
            return i
            #I.append(i)
    #return I
    return -1
    

def intersect_samples(samples1, samples2, case_insensitive):
    # Return the intersection of samples1 and samples2.  Preserves the
    # order according to samples1.

    ## TOO SLOW.
    ## in_both = []
    ## for x in annots1:
    ##     matching = False
    ##     for y in annots2:
    ##         if is_same_annot(x, y, case_insensitive):
    ##             matching = True
    ##             break
    ##     if matching:
    ##         in_both.append(x)
    ## return in_both
    
    samples1_cmp = samples1
    samples2_cmp = samples2
    if case_insensitive:
        samples1_cmp = [x.upper() for x in samples1]
        samples2_cmp = [x.upper() for x in samples2]
    samples2_cmp = {}.fromkeys(samples2_cmp)
    in_both = [
        samples1[i] for i in range(len(samples1)) if
        samples1_cmp[i] in samples2_cmp]
    return in_both


def write_matrix(outfile, matrix):
    if isinstance(matrix, AnnotationMatrix):
        write_annot(outfile, matrix)
    else:
        write_express(outfile, matrix)


def read_express(filename):
    import arrayio

    # Read the files.
    assert os.path.exists(filename), \
        "I could not find the file: %s" % filename
    fmt_module = arrayio.choose_format(filename)
    assert fmt_module, \
        "I could not figure out the format of file: %s" % filename
    x = fmt_module.read(filename)
    return x


def write_express(filename, matrix):
    import arrayio
    arrayio.write(matrix, open(filename, 'w'))


def read_annot(filename):
    from genomicode import genesetlib

    name_order = []
    name2annots = {}
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, annots = x
        name_order.append(name)
        name2annots[name] = annots
    return AnnotationMatrix(name2annots, name_order)


def write_annot(filename, annot_matrix):
    from genomicode import jmath
    matrix = []
    for name in annot_matrix.name_order:
        annots = annot_matrix.name2annots[name]
        x = [name] + annots
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)

    handle = open(filename, 'w')
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Align a set of matrices.  Preserve the order of the "
        "first file given.")
    parser.add_argument("outfile", nargs="+")
    
    parser.add_argument(
        "--clobber", default=False, action="store_true",
        help="Overwrite output files, if they already exist.")
    parser.add_argument(
        "--case_insensitive", default=False, action="store_true",
        help="Do a case insensitive search of sample names.")
    parser.add_argument(
        "--outer_join", default=False, action="store_true",
        help='By default, does an "inner join" and keeps only the '
        'records that are present in all files.  An "outer join" will '
        'also keep records that occur in one file.')
    parser.add_argument(
        "--null_string", default="",
        help='For outer_join, what to give the missing values.')

    parser.add_argument(
        "--express_file", default=[], action="append", help="")
    parser.add_argument(
        "--annot_file", default=[], action="append", help="")

    #parser.add_argument(
    #    "--headers", help="Explicitly specify the headers for each of the "
    #    "file.  The headers should be given as a semicolon-separated list.")
    parser.add_argument(
        "--first_annot_header", help="If only aligning annotation files, "
        "find the samples to be matched under this header in the first "
        "annotation file.")
    
    args = parser.parse_args()
    assert len(args.outfile) == len(args.express_file) + len(args.annot_file)
    for x in args.express_file + args.annot_file:
        assert os.path.exists(x), "I could not find file: %s" % x
    for x in args.outfile:
        assert args.clobber or not os.path.exists(x), "File exists: %s" % x
    if args.null_string:
        assert args.outer_join, \
               "null_string given, but only used for outer_join"

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
    new_matrix_data = []  # list of (infile, outfile, matrix)
    for x in matrix_data:
        infile, outfile, is_express_file = x
        if is_express_file:
            data = read_express(infile)
        else:
            data = read_annot(infile)
        x = infile, outfile, data
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    if args.outer_join:
        samples = list_all_samples(
            matrix_data, args.first_annot_header, args.case_insensitive)
        assert samples, "No samples."
    else:
        samples = list_common_samples(
            matrix_data, args.first_annot_header, args.case_insensitive)
        assert samples, "No common samples found."

    # Align each of the matrices.
    matrix_data = align_matrices(
        matrix_data, samples, args.case_insensitive, args.null_string)

    # Write out each of the matrices.
    for x in matrix_data:
        infile, outfile, matrix = x
        write_matrix(outfile, matrix)
    

if __name__ == '__main__':
    main()
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
