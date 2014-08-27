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
# uniq_samples
# strip_nonalnum
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


def list_all_samples(matrix_data, case_insensitive, hash_samples,
                     ignore_nonalnum):
    from genomicode import hashlib
    
    assert matrix_data

    # Get the samples that occur in any of the files.  Preserve the
    # order of the samples.
    all_samples = []
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x
        all_samples.extend(samples)

    # Get rid of duplicate samples.
    duplicated = [False] * len(all_samples)
    all_samples_cmp = all_samples
    if case_insensitive:
        all_samples_cmp = [x.upper() for x in all_samples_cmp]
    if hash_samples:
        all_samples_cmp = [hashlib.hash_var(x) for x in all_samples_cmp]
    if ignore_nonalnum:
        all_samples_cmp = [strip_nonalnum(x) for x in all_samples_cmp]

    seen = {}
    for i, s in enumerate(all_samples_cmp):
        if s in seen:
            duplicated[i] = True
        seen[s] = 1

    all_samples = [
        all_samples[i] for i in range(len(all_samples)) if not duplicated[i]]
    return all_samples


def list_common_samples(matrix_data, case_insensitive, hash_samples,
                        ignore_nonalnum):
    assert matrix_data

    # Initialize the common samples.
    samples_hint = peek_samples_hint(matrix_data)
    assert samples_hint
    common = samples_hint

    # Get the samples that are in all of the data sets.  Keep track of
    # the order of the samples in the first file so we can resort them
    # later.
    first_samples = None
    for i, x in enumerate(matrix_data):
        infile, outfile, matrix, header, samples = x
        common = intersect_samples(
            common, samples, case_insensitive, hash_samples, ignore_nonalnum)
        if i == 0:
            first_samples = samples

    # Finally, order the annotations according to the first file given.  
    assert first_samples
    samples = intersect_samples(
        first_samples, common, case_insensitive, hash_samples, ignore_nonalnum)
    
    return samples


def peek_samples_hint(matrix_data):
    # Figure out what the samples look like.  Look for an expression
    # file first.  If an expression file is not found, use the
    # annotation files.
    samples_hint = None
    x = [x for x in matrix_data if not isinstance(x[2], AnnotationMatrix)]
    if x:
        # Found an expression file.
        infile, outfile, matrix, header = x[0][:4]
        assert header is None
        samples_hint = get_express_samples(matrix)
    else:
        # No expression files.  Look for an annotation file with given
        # --header.
        samples_hint = None
        for x in matrix_data:
            infile, outfile, matrix, header = matrix_data[0][:4]
            if header is None:
                continue
            assert header in matrix.name2annots, "Missing header: %s\n%s" % (
                repr(header), sorted(matrix.name2annots))
            samples_hint = matrix.name2annots[header]
            break
        assert samples_hint is not None, \
               "Please specify at least one --header."
    assert samples_hint
    return samples_hint


def align_matrices(
    matrix_data, align_samples, case_insensitive, hash_samples,
    ignore_nonalnum, ignore_blank, left_join, outer_join, null_string):
    # align_samples is a list of unique samples.
    import itertools
    from genomicode import hashlib

    sample2matrix2indexes = {}  # sample_i -> matrix_i -> list of indexes
    for i, sample in enumerate(align_samples):
        sample2matrix2indexes[i] = {}
        for j, x in enumerate(matrix_data):
            infile, outfile, matrix, header, samples = x

            samples_cmp = samples
            sample_cmp = sample
            if case_insensitive:
                samples_cmp = [x.upper() for x in samples_cmp]
                sample_cmp = sample_cmp.upper()
            if hash_samples:
                samples_cmp = [hashlib.hash_var(x) for x in samples_cmp]
                sample_cmp = hashlib.hash_var(sample_cmp)
            if ignore_nonalnum:
                samples_cmp = [strip_nonalnum(x) for x in samples_cmp]
                sample_cmp = strip_nonalnum(sample_cmp)

            # Keep the blanks in the first matrix.  Just don't align
            # them to anything in the remaining matrices.
            if ignore_blank and not sample_cmp.strip() and j > 0:
                indexes = []
            else:
                indexes = [
                    k for k in range(len(samples_cmp))
                    if samples_cmp[k] == sample_cmp]
            sample2matrix2indexes[i][j] = indexes

    # Now align the indexes for each matrix.  Here, len(indexes)
    # should be the same for each matrix.
    sample2matrix2aligned = {}  # sample_i -> matrix_i -> list of indexes
    for i in range(len(align_samples)):
        matrix2indexes = sample2matrix2indexes[i].copy()
        matrix2aligned = {}

        if left_join or outer_join:
            for j in range(len(matrix_data)):
                matrix2aligned[j] = []

            # list (of len(matrix_data)) of list of indexes
            all_indexes = []  
            for j in range(len(matrix_data)):
                indexes = matrix2indexes[j]
                if not indexes:
                    indexes = [None]
                all_indexes.append(indexes)

            # Do a product of all matching indexes across the matrices.
            for x in itertools.product(*all_indexes):
                assert len(x) == len(matrix_data)
                for j in range(len(matrix_data)):
                    matrix2aligned[j].append(x[j])
        else:
            # If this sample is not found in all data sets, then
            # ignore it.
            found = True
            for j in range(len(matrix_data)):
                if not matrix2indexes[j]:
                    found = False
            if not found:
                for j in range(len(matrix_data)):
                    matrix2aligned[j] = []
            else:
                # Use only the first instance of each sample.
                for j in range(len(matrix_data)):
                    indexes = matrix2indexes[j]
                    matrix2aligned[j] = [indexes[0]]
            
        sample2matrix2aligned[i] = matrix2aligned

    # Make the list of indexes for each matrix.
    matrix2indexes = {}
    for j in range(len(matrix_data)):
        matrix2indexes[j] = []
    for i in range(len(align_samples)):
        for j in range(len(matrix_data)):
            indexes = sample2matrix2aligned[i][j]
            matrix2indexes[j].extend(indexes)

    #for i in range(len(matrix2indexes[0])):
    #    x = [matrix2indexes[j][i] for j in range(len(matrix_data))]
    #    print "\t".join(map(str, x))
    #print [len(matrix2indexes[i]) for i in range(len(matrix_data))]
    #import sys; sys.exit(0)

    aligned_matrix_data = []
    for j, x in enumerate(matrix_data):
        infile, outfile, matrix, header, samples = x
        aligned_matrix = align_matrix(matrix, matrix2indexes[j], null_string)
        x = infile, outfile, aligned_matrix, header, samples
        aligned_matrix_data.append(x)
    return aligned_matrix_data
            
            
def align_matrix(matrix, indexes, null_string):
    if isinstance(matrix, AnnotationMatrix):
        aligned_matrix = align_annot(matrix, indexes, null_string)
    else:
        aligned_matrix = align_express(matrix, indexes, null_string)
    return aligned_matrix


def align_express(matrix, indexes, null_string):
    import arrayio

    # The Matrix class requires actual indexes and not None.  Give it
    # index 0, and then remove those lines myself later.
    I_index = []
    for i in indexes:
        if i != None:
            I_index.append(i)
        else:
            I_index.append(0)
    matrix_aligned = matrix.matrix(None, I_index)

    # Fix the None indexes in the matrix.
    X = matrix_aligned._X
    assert len(indexes) == len(X[0])
    for i_new, i_old in enumerate(indexes):
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
        for i_new, i_old in enumerate(indexes):
            if i_old != None:
                annots_new.append(annots[i_new])
            elif name == header:
                # Should figure out how to pull out the new sample
                # names.
                #annots_new.append(samples[i_new])
                annots_new.append(null_string)
            else:
                annots_new.append(null_string)
        matrix_aligned._col_names[name] = annots_new
    
    return matrix_aligned


def align_annot(matrix, indexes, null_string):
    name2annots_new = {}
    for name, annots in matrix.name2annots.iteritems():
        annots_new = []
        for i, i_annot in enumerate(indexes):
            if i_annot != None:
                annots_new.append(annots[i_annot])
            #elif name == header:
            #    annots_new.append(samples[i])
            else:
                annots_new.append(null_string)
        name2annots_new[name] = annots_new
    return AnnotationMatrix(name2annots_new, matrix.name_order)


## def find_matrix_samples(matrix, header, sample, case_insensitive):
##     if isinstance(matrix, AnnotationMatrix):
##         return find_annot_samples(matrix, header, sample, case_insensitive)
##     return find_express_samples(matrix, header, sample, case_insensitive)


## def find_annot_samples(matrix, header, sample, case_insensitive):
##     assert header in matrix.name2annots, header

##     sample_list = matrix.name2annots[header]

##     sample_list_cmp = sample_list
##     sample_cmp = sample
##     if case_insensitive:
##         sample_list_cmp = [x.upper() for x in sample_list]
##         sample_cmp = sample.upper()

##     indexes = [
##         i for i in range(len(sample_list_cmp))
##         if sample_list_cmp[i] == sample_cmp]
##     return indexes

    
## def find_express_samples(matrix, header, sample, case_insensitive):
##     import arrayio
    
##     name = arrayio.COL_ID
##     if name not in matrix._col_names:
##         name = matrix._synonyms[name]
##     assert name in matrix._col_names, "I can not find the sample names."
##     sample_list = matrix.col_names(name)

##     sample_list_cmp = sample_list
##     sample_cmp = sample
##     if case_insensitive:
##         sample_list_cmp = [x.upper() for x in sample_list]
##         sample_cmp = sample.upper()

##     indexes = [
##         i for i in range(len(sample_list_cmp))
##         if sample_list_cmp[i] == sample_cmp]
##     return indexes

    
## def align_matrices(matrix_data, samples, case_insensitive, null_string):
##     new_matrix_data = []
##     for x in matrix_data:
##         infile, outfile, matrix, header = x
##         aligned_matrix = align_matrix(
##             matrix, header, samples, case_insensitive, null_string)
##         x = infile, outfile, aligned_matrix, header
##         new_matrix_data.append(x)
##     return new_matrix_data


## def align_matrix(
##     matrix, header, samples, case_insensitive, null_string):
##     if isinstance(matrix, AnnotationMatrix):
##         aligned_matrix = align_annot(
##             matrix, header, samples, case_insensitive, null_string)
##     else:
##         assert header is None
##         aligned_matrix = align_express(
##             matrix, samples, case_insensitive, null_string)
##     return aligned_matrix


## def align_express(matrix, samples, case_insensitive, null_string):
##     import arrayio
    
##     names = get_express_samples(matrix)
##     I = []
##     for n in samples:
##         i = find_sample(names, n, case_insensitive)
##         if i == -1:
##             i = None
##         I.append(i)

##     # The Matrix class requires actually indexes.  Give it index 0,
##     # and then remove those lines myself later.
##     I_index = []
##     for i in I:
##         if i != None:
##             I_index.append(i)
##         else:
##             I_index.append(0)
##     assert len(I_index) == len(I)
##     matrix_aligned = matrix.matrix(None, I_index)

##     # Fix the None indexes in the matrix.
##     X = matrix_aligned._X
##     assert len(I) == len(X[0])
##     for i_new, i_old in enumerate(I):
##         if i_old != None:
##             continue
##         # Fill column i with blanks.
##         for j in range(len(X)):
##             X[j][i_new] = ""

##     # Fix the annotations.
##     header = arrayio.COL_ID
##     if header not in matrix._col_names:
##         header = matrix._synonyms[header]
##     assert header in matrix._col_names, "I can not find the sample names."
    
##     # Fix the sample names.
##     for name in matrix_aligned._col_names:
##         annots = matrix_aligned._col_names[name]
##         annots_new = []
##         for i_new, i_old in enumerate(I):
##             if i_old != None:
##                 annots_new.append(annots[i_new])
##             elif name == header:
##                 annots_new.append(samples[i_new])
##             else:
##                 annots_new.append(null_string)
##         matrix_aligned._col_names[name] = annots_new
    
##     return matrix_aligned


## def align_annot(matrix, header, samples, case_insensitive, null_string):
##     names = get_annot_samples(matrix, header, samples, case_insensitive)

##     # Figure out which is the header for these samples.
##     header = None
##     for h, n in matrix.name2annots.iteritems():
##         if n == names:
##             header = h
##     assert header is not None

##     # Index each name in the matrix.
##     ## Too slow.
##     ##I = []
##     ##for n in samples_cmp:
##     ##    i = find_sample(names, n, case_insensitive)
##     ##    if i == -1:
##     ##        i = None
##     ##    I.append(i)
##     ##assert len(I) == len(samples)
##     names_cmp = names
##     samples_cmp = samples
##     if case_insensitive:
##         names_cmp = [x.upper() for x in names]
##         samples_cmp = [x.upper() for x in samples]
##     name2indexes = {}  # name -> list of indexes
##     for (i, n) in enumerate(names_cmp):
##         if n not in name2indexes:
##             name2indexes[n] = []
##         name2indexes[n].append(i)

##     # Find the indexes of each sample.  If a sample is requested
##     # multiple times, and also occurs multiple times in the file,
##     # return distinct records from the file.
##     name2inum = {}
##     for n in name2indexes:
##         name2inum[n] = 0
##     I = []
##     for n in samples_cmp:
##         if n not in name2indexes:
##             I.append(None)
##             continue
##         indexes = name2indexes[n]
##         i = name2inum[n]
##         if i >= len(indexes):
##             i = 0
##         name2inum[n] = i+1
##         I.append(indexes[i])
        
##     name2annots_new = {}
##     for name, annots in matrix.name2annots.iteritems():
##         annots_new = []
##         for i, i_annot in enumerate(I):
##             if i_annot != None:
##                 annots_new.append(annots[i_annot])
##             elif name == header:
##                 annots_new.append(samples[i])
##             else:
##                 annots_new.append(null_string)
##         name2annots_new[name] = annots_new
##     return AnnotationMatrix(name2annots_new, matrix.name_order)


def get_samples(
    matrix, header_hint, samples_hint, case_insensitive, hash_samples,
    ignore_nonalnum):
    # Get the samples from the matrix.  Since in principle anything in
    # the annotation file can be a sample, need to give it a hint of
    # what the samples look like.
    if isinstance(matrix, AnnotationMatrix):
        return get_annot_samples(
            matrix, header_hint, samples_hint, case_insensitive, hash_samples,
            ignore_nonalnum)
    return get_express_samples(matrix)


def get_express_samples(matrix):
    import arrayio
    
    name = arrayio.COL_ID
    if name not in matrix._col_names:
        name = matrix._synonyms[name]
    assert name in matrix._col_names, "I can not find the sample names."
    x = matrix.col_names(name)
    return x


def get_annot_samples(matrix, header_hint, samples_hint,
                      case_insensitive, hash_samples, ignore_nonalnum):
    if header_hint:
        assert header_hint in matrix.name2annots
        return matrix.name2annots[header_hint]
    
    all_matches = []  # list of (num_matches, name, matches)
    for name, annots in matrix.name2annots.iteritems():
        x = intersect_samples(
            annots, samples_hint, case_insensitive, hash_samples,
            ignore_nonalnum)
        matches = uniq_samples(
            x, case_insensitive, hash_samples, ignore_nonalnum)
        x = len(matches), name, matches
        all_matches.append(x)
    assert all_matches
    all_matches = sorted(all_matches)
    x = all_matches[-1]
    x, name, x = x
    return matrix.name2annots[name]


def cmp_sample(x, y,
               case_insensitive, hash_samples, ignore_nonalnum, ignore_blank):
    from genomicode import hashlib

    if ignore_blank and x.strip():
        return False
    
    if case_insensitive:
        x, y = x.upper(), y.upper()
    if hash_samples:
        x, y = hashlib.hash_var(x), hashlib.hash_var(y)
    if ignore_nonalnum:
        x, y = strip_nonalnum(x), strip_nonalnum(y)
    return x == y


## def index_sample(sample_list, sample, case_insensitive):
##     i = find_sample(sample_list, sample, case_insensitive)
##     assert i >= 0, "Missing: %s" % sample


def find_sample(sample_list, sample, case_insensitive, hash_samples,
                ignore_nonalnum, ignore_blank):
    # Return the index of this sample or -1.
    from genomicode import hashlib

    if ignore_blank and not sample.strip():
        return -1

    sample_list_cmp = sample_list
    sample_cmp = sample
    if case_insensitive:
        sample_list_cmp = [x.upper() for x in sample_list_cmp]
        sample_cmp = sample_cmp.upper()
    if hash_samples:
        sample_list_cmp = [hashlib.hash_var(x) for x in sample_list_cmp]
        sample_cmp = hashlib.hash_var(sample_cmp)
    if ignore_nonalnum:
        sample_list_cmp = [strip_nonalnum(x) for x in sample_list_cmp]
        sample_cmp = strip_nonalnum(sample_cmp)

    ## Too slow.
    ##for i, x in enumerate(sample_list_cmp):
    ##    if x == sample_cmp:
    ##        return i
    ##        #I.append(i)
    ##return -1
    try:
        return sample_list_cmp.index(sample_cmp)
    except ValueError, x:
        pass
    return -1
    

def intersect_samples(samples1, samples2, case_insensitive, hash_samples,
                      ignore_nonalnum):
    # Return the intersection of samples1 and samples2.  Preserves the
    # order according to samples1.
    from genomicode import hashlib

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
        samples1_cmp = [x.upper() for x in samples1_cmp]
        samples2_cmp = [x.upper() for x in samples2_cmp]
    if hash_samples:
        samples1_cmp = [hashlib.hash_var(x) for x in samples1_cmp]
        samples2_cmp = [hashlib.hash_var(x) for x in samples2_cmp]
    if ignore_nonalnum:
        samples1_cmp = [strip_nonalnum(x) for x in samples1_cmp]
        samples2_cmp = [strip_nonalnum(x) for x in samples2_cmp]

    samples2_cmp = {}.fromkeys(samples2_cmp)
    in_both = [
        samples1[i] for i in range(len(samples1)) if
        samples1_cmp[i] in samples2_cmp]
    return in_both


def uniq_samples(samples, case_insensitive, hash_samples, ignore_nonalnum):
    from genomicode import hashlib

    samples_cmp = samples
    if case_insensitive:
        samples_cmp = [x.upper() for x in samples_cmp]
    if hash_samples:
        samples_cmp = [hashlib.hash_var(x) for x in samples_cmp]
    if ignore_nonalnum:
        samples_cmp = [strip_nonalnum(x) for x in samples_cmp]
        
    seen = {}  # dict of sample_cmp -> sample
    for sample, sample_cmp in zip(samples, samples_cmp):
        if sample_cmp in seen:
            continue
        seen[sample_cmp] = sample
        
    # Return a list of the unique samples.
    return seen.values()


def strip_nonalnum(s):
    import re
    s = re.sub("\W", "", s)
    s = s.replace("_", "")
    return s


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
    import re
    from genomicode import genesetlib

    name_order = []
    name2annots = {}
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, annots = x

        # Hack: Some files contain special characters, which mess up
        # alignment. Fix this here.
        # na\xc3\xafve-WIBR3.5 hESC
        # na\xe2\x80\x9a\xc3\xa0\xc3\xb6\xe2\x88\x9a\xc3\xb2ve-C1.2 hiPSC
        annots = [re.sub("na\\W+ve", "naive", x) for x in annots]
        
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
        "--strict", default=False, action="store_true",
        help="Complain if a file is missing a sample.")
    parser.add_argument(
        "--case_insensitive", default=False, action="store_true",
        help="Do a case insensitive search of sample names.")
    parser.add_argument(
        "--hash", default=False, action="store_true",
        help="Hash the sample names to [a-zA-Z0-9_] before comparison.")
    parser.add_argument(
        "--ignore_nonalnum", default=False, action="store_true",
        help="Ignore non-alphanumeric characters in the IDs.")
    parser.add_argument(
        "--ignore_blank", default=False, action="store_true",
        help="Ignore IDs that are blank (don't align them.")
    parser.add_argument(
        "--left_join", default=False, action="store_true",
        help='By default, does an "inner join" and keeps only the '
        'records that are present in all files.  A "left join" will '
        'keep all records that occur in the first file.')
    parser.add_argument(
        "--outer_join", default=False, action="store_true",
        help='By default, does an "inner join" and keeps only the '
        'records that are present in all files.  An "outer join" will '
        'also keep records that occur in any file.')
    parser.add_argument(
        "--null_string", default="",
        help='For left_join or outer_join, what to give the missing values.')

    parser.add_argument(
        "--express_file", default=[], action="append", help="")
    parser.add_argument(
        "--annot_file", default=[], action="append", help="")

    parser.add_argument(
        "--header", default=[], action="append",
        help="Specify the header for an annotation file.")
    #parser.add_argument(
    #    "--first_annot_header", help="If only aligning annotation files, "
    #    "find the samples to be matched under this header in the first "
    #    "annotation file.")
    
    args = parser.parse_args()
    assert len(args.outfile) == len(args.express_file) + len(args.annot_file)
    for x in args.express_file + args.annot_file:
        assert os.path.exists(x), "I could not find file: %s" % x
    for x in args.outfile:
        assert args.clobber or not os.path.exists(x), "File exists: %s" % x
    assert not (args.left_join and args.outer_join)
    if args.null_string:
        assert args.outer_join or args.left_join, \
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

    # Align the headers to the annotation files.
    headers = [None] * len(matrix_data)
    header_i = -1
    for i, arg in enumerate(sys.argv):
        if arg == "--header":
            assert header_i >= 0, \
                   "--header given before an --express_file or --annot_file."
            assert headers[header_i] is None, "Two --header for one file."
            headers[header_i] = sys.argv[i+1]
        elif arg in ["--express_file", "--annot_file"]:
            header_i += 1

    # Add the headers to the matrix_data.
    new_matrix_data = []  # list of (infile, outfile, is_express_file, header)
    for i in range(len(matrix_data)):
        infile, outfile, is_express_file = matrix_data[i]
        if is_express_file and headers[i]:
            raise NotImplementedError, "No headers for --express_file."
        x = infile, outfile, is_express_file, headers[i]
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    # Read each of the files.
    new_matrix_data = []  # list of (infile, outfile, matrix, header)
    for x in matrix_data:
        infile, outfile, is_express_file, header = x
        if is_express_file:
            data = read_express(infile)
        else:
            data = read_annot(infile)
        x = infile, outfile, data, header
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    # Find the samples in each matrix.
    new_matrix_data = []  # list of (infile, outfile, matrix, header, samples)
    samples_hint = peek_samples_hint(matrix_data)
    for x in matrix_data:
        infile, outfile, matrix, header = x
        samples = get_samples(
            matrix, header, samples_hint, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        x = infile, outfile, matrix, header, samples
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    if args.left_join:
        assert not args.strict, "Can't do a strict left join."
        samples = list_all_samples(
            matrix_data[:1], args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        assert samples, "No samples."
    elif args.outer_join:
        assert not args.strict, "Can't do a strict outer join."
        samples = list_all_samples(
            matrix_data, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        assert samples, "No samples."
    else:
        samples = list_common_samples(
            matrix_data, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        assert samples, "No common samples found."

    if args.strict:
        all_samples = list_all_samples(
            matrix_data, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        common_samples = list_common_samples(
            matrix_data, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        if sorted(all_samples) != sorted(common_samples):
            missing_samples = []
            for x in all_samples:
                i = find_sample(
                    common_samples, x, args.case_insensitive, args.hash,
                    args.ignore_nonalnum, args.ignore_blank)
                if i >= 0:
                    continue
                missing_samples.append(x)
            short = missing_samples
            if len(short) > 10:
                short = short[:10] + ["..."]
            short = "\n".join(short)
            raise AssertionError, "%d samples not in all data sets.\n%s" % \
                  (len(missing_samples), short)

    # Align each of the matrices.
    matrix_data = align_matrices(
        matrix_data, samples, args.case_insensitive, args.hash,
        args.ignore_nonalnum, args.ignore_blank,
        args.left_join, args.outer_join, args.null_string)

    # Write out each of the matrices.
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x
        write_matrix(outfile, matrix)
    

if __name__ == '__main__':
    main()
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
