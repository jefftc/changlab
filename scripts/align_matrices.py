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
# align_samples
#
# add_missing_samples
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
        #assert header is None
        x, samples_hint = get_express_samples(matrix)
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


def _process_sample(sample, case_insensitive, hash_samples, ignore_nonalnum):
    from genomicode import hashlib

    x = sample
    if case_insensitive:
        x = x.upper()
    if hash_samples:
        x = hashlib.hash_var(x)
    if ignore_nonalnum:
        x = strip_nonalnum(x)
    return x


def _left_join_matrices(matrix_samples_cmp, all_samples,
                        sample2matrix2indexes):
    # Return matrix2indexes as a list of lists.  The outer list
    # corresponds to each matrix.  The inner is a list of indexes
    # that indicate the index into the matrix.  Any value can be None,
    # which indicates the value is not in the matrix.

    # all_samples               list of samples
    # matrix_samples_cmp        list of samples
    # sample2matrix2indexes     sample_i -> matrix_i -> list of indexes
    assert matrix_samples_cmp
    matrix2indexes = [[] for i in range(len(matrix_samples_cmp))]
    sample2n = {}  # sample_i -> number of times this sample already seen.

    matrix_samples = matrix_samples_cmp[0]
    for sample in matrix_samples:
        sample_i = all_samples.index(sample)
        if sample_i not in sample2n:
            sample2n[sample_i] = 0
        for i in range(len(matrix_samples_cmp)):
            indexes = sample2matrix2indexes[sample_i][i]
            # For left_join, take all the samples from the first
            # matrix, and the first matching sample from the
            # subsequent ones.
            if i == 0:
                assert sample2n[sample_i] < len(indexes)
                index = indexes[sample2n[sample_i]]
            elif indexes:
                index = indexes[0]
            else:
                index = None
            matrix2indexes[i].append(index)
        sample2n[sample_i] += 1

    return matrix2indexes


def _outer_join_matrices(matrix_samples_cmp, all_samples,
                         sample2matrix2indexes):
    # Return matrix2indexes as a list of lists.  The first list
    # corresponds to each matrix.  The second is a list of indexes
    # that indicate the index into the matrix.  Any value can be None,
    # which indicates the value is not in the matrix.
    import itertools

    # all_samples               list of samples
    # matrix_samples_cmp        list of list of samples
    # sample2matrix2indexes     sample_i -> matrix_i -> list of indexes
    assert matrix_samples_cmp
    matrix2indexes = [[] for i in range(len(matrix_samples_cmp))]

    for sample_i in range(len(all_samples)):
        # Do a product of all matching indexes across the matrices.
        all_indexes = []
        for i in range(len(matrix_samples_cmp)):
            x = sample2matrix2indexes[sample_i][i]
            if not x:
                x = [None]
            all_indexes.append(x)
        for indexes in itertools.product(*all_indexes):
            assert len(indexes) == len(matrix_samples_cmp)
            for i in range(len(indexes)):
                matrix2indexes[i].append(indexes[i])
    return matrix2indexes


def _inner_join_matrices(matrix_samples_cmp, all_samples,
                         sample2matrix2indexes):
    # Return matrix2indexes as a list of lists.  The first list
    # corresponds to each matrix.  The second is a list of indexes
    # that indicate the index into the matrix.
    import itertools

    # all_samples               list of samples
    # matrix_samples_cmp        list of list of samples
    # sample2matrix2indexes     sample_i -> matrix_i -> list of indexes
    assert matrix_samples_cmp
    matrix2indexes = [[] for i in range(len(matrix_samples_cmp))]

    for sample_i in range(len(all_samples)):
        # Take the first sample that matches.
        for i in range(len(matrix_samples_cmp)):
            x = sample2matrix2indexes[sample_i][i]
            assert x
            matrix2indexes[i].append(x[0])
    return matrix2indexes


def align_matrices(
    matrix_data, final_samples, case_insensitive, hash_samples,
    ignore_nonalnum, ignore_blank, left_join, outer_join, unaligned_only,
    null_string):
    # final_samples is a list of unique samples to be included in the
    # final matrix.
    import itertools
    from genomicode import hashlib

    final_samples_cmp = [
        _process_sample(x, case_insensitive, hash_samples, ignore_nonalnum)
        for x in final_samples]

    # Pre-process the samples so I don't have to do it repeatedly.
    # This list should be aligned to matrix_data.
    matrix_samples_cmp = []
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x
        x = [
            _process_sample(x, case_insensitive, hash_samples, ignore_nonalnum)
            for x in samples]
        matrix_samples_cmp.append(x)

    # Pre-process, for each matrix, a dictionary of sample -> list of
    # indexes.
    matrix_sample2indexes = []
    for samples_cmp in matrix_samples_cmp:
        sample2indexes = {}
        for i, s in enumerate(samples_cmp):
            if s not in sample2indexes:
                sample2indexes[s] = []
            sample2indexes[s].append(i)
        matrix_sample2indexes.append(sample2indexes)

    # Reorganize data structure.  Want to be able to access location
    # of each sample.
    sample2matrix2indexes = {}  # sample_i -> matrix_i -> list of indexes
    for i, sample in enumerate(final_samples):
        sample2matrix2indexes[i] = {}
        sample_cmp = final_samples_cmp[i]
        for j, x in enumerate(matrix_data):
            sample2indexes = matrix_sample2indexes[j]

            # Keep the blanks in the first matrix.  Just don't align
            # them to anything in the remaining matrices.
            if ignore_blank and not sample_cmp.strip() and j > 0:
                indexes = []
            else:
                indexes = sample2indexes.get(sample_cmp, [])
            sample2matrix2indexes[i][j] = indexes

    # Now align the indexes for each matrix.  Here, len(indexes)
    # should be the same for each matrix.
    if left_join:
        matrix2indexes = _left_join_matrices(
            matrix_samples_cmp, final_samples, sample2matrix2indexes)
    elif outer_join:
        matrix2indexes = _outer_join_matrices(
            matrix_samples_cmp, final_samples, sample2matrix2indexes)
    else:
        matrix2indexes = _inner_join_matrices(
            matrix_samples_cmp, final_samples, sample2matrix2indexes)

    if unaligned_only:
        # Make a list of all the indexes that are unaligned.
        I = []  # unaligned indexes
        for indexes in matrix2indexes:
            for i, index in enumerate(indexes):
                if index is None:
                    I.append(i)
        I = sorted({}.fromkeys(I))
        # Pull out just the indexes that are unaligned.
        for i in range(len(matrix2indexes)):
            x = matrix2indexes[i]
            x = [x[j] for j in I]
            matrix2indexes[i] = x
        

    aligned_matrix_data = []
    for j, x in enumerate(matrix_data):
        infile, outfile, matrix, header, samples = x
        aligned_matrix = align_matrix(matrix, matrix2indexes[j], null_string)
        aligned_samples = align_samples(
            samples, matrix2indexes[j], null_string)
        x = infile, outfile, aligned_matrix, header, aligned_samples
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

    # The Matrix class requires actual indexes and not None.  For
    # None, give it index 0, and then remove those lines myself later.
    I_index = []
    for i in indexes:
        if i != None:
            I_index.append(i)
        else:
            I_index.append(0)
    matrix_aligned = matrix.matrix(None, I_index)

    # Fix the None indexes in the matrix.
    X = matrix_aligned._X
    if indexes:
        assert len(indexes) == len(X[0])
    else:
        assert not X
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


def align_samples(samples, indexes, null_string):
    samples_new = []
    for i, i_annot in enumerate(indexes):
        if i_annot != None:
            samples_new.append(samples[i_annot])
        else:
            samples_new.append(null_string)
    return samples_new


def add_missing_samples(matrix_data, null_string):
    # Make a list of all the samples.
    complete_samples = None
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x
        if not complete_samples:
            complete_samples = samples
        assert len(complete_samples) == len(samples)
        for i in range(len(complete_samples)):
            if complete_samples[i] == null_string:
                complete_samples[i] = samples[i]

    # Now fix each of the matrices.
    new_matrix_data = []
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x

        if isinstance(matrix, AnnotationMatrix):
            samples = matrix.name2annots[header]
        else:
            samples = matrix.col_names(header)
        for i in range(len(samples)):
            if samples[i] == null_string:
                samples[i] = complete_samples[i]
        if isinstance(matrix, AnnotationMatrix):
            matrix.name2annots[header] = samples
        else:
            matrix._col_names[header] = samples
        x = infile, outfile, matrix, header, samples
        new_matrix_data.append(x)
        
    return new_matrix_data


def get_samples(
    matrix, header_hint, samples_hint, case_insensitive, hash_samples,
    ignore_nonalnum):
    # Return (header, samples) from the matrix.  Since in principle
    # anything in the annotation file can be a sample, need to give it
    # a hint of what the samples look like.  Returns None if I could
    # not find the samples.
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
    return name, x


def get_annot_samples(matrix, header_hint, samples_hint,
                      case_insensitive, hash_samples, ignore_nonalnum):
    if header_hint:
        assert header_hint in matrix.name2annots
        return header_hint, matrix.name2annots[header_hint]

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
    num_matches, name, x = x
    if not num_matches:
        return None
    return name, matrix.name2annots[name]


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

    SKIP_OUTFILE = "_"

    parser = argparse.ArgumentParser(
        description="Align a set of matrices.  Preserve the order of the "
        "first file given.")
    parser.add_argument("outfile", nargs="+")

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
    parser.add_argument(
        "--clobber", default=False, action="store_true",
        help="Overwrite output files, if they already exist.")

    group = parser.add_argument_group(title="Comparisons")
    group.add_argument(
        "--case_insensitive", default=False, action="store_true",
        help="Do a case insensitive search of sample names.")
    group.add_argument(
        "--hash", default=False, action="store_true",
        help="Hash the sample names to [a-zA-Z0-9_] before comparison.")
    group.add_argument(
        "--ignore_nonalnum", default=False, action="store_true",
        help="Ignore non-alphanumeric characters in the IDs.")
    group.add_argument(
        "--ignore_blank", default=False, action="store_true",
        help="Ignore IDs that are blank (don't align them.")

    group = parser.add_argument_group(title="Joins")
    group.add_argument(
        "--strict", default=False, action="store_true",
        help="Complain if a file is missing a sample.")
    group.add_argument(
        "--left_join", default=False, action="store_true",
        help='By default, does an "inner join" and keeps only the '
        'records that are present in all files.  A "left join" will '
        'keep all records that occur in the first file.')
    group.add_argument(
        "--outer_join", default=False, action="store_true",
        help='By default, does an "inner join" and keeps only the '
        'records that are present in all files.  An "outer join" will '
        'also keep records that occur in any file.')


    group = parser.add_argument_group(title="Output")
    group.add_argument(
        "--null_string", default="",
        help='For left_join or outer_join, what to give the missing values.')
    group.add_argument(
        "--unaligned_only", action="store_true",
        help="Show only the rows that are not aligned.")
    group.add_argument(
        "--dont_add_missing_samples", action="store_true",
        help="If a matrix does not have a sample, don't fill in the value "
        "from another matrix.")


    args = parser.parse_args()
    ni, no = len(args.express_file)+len(args.annot_file), len(args.outfile)
    assert ni == no, "Mismatch: %d inputs and %d outputs" % (ni, no)
        
    for x in args.express_file + args.annot_file:
        assert os.path.exists(x), "I could not find file: %s" % x
    for x in args.outfile:
        if x == SKIP_OUTFILE:
            continue
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
        x = get_samples(
            matrix, header, samples_hint, args.case_insensitive, args.hash,
            args.ignore_nonalnum)
        assert x, "I could not find the samples for %s" % infile
        header, samples = x
        x = infile, outfile, matrix, header, samples
        new_matrix_data.append(x)
    matrix_data = new_matrix_data

    if args.left_join:
        assert not args.strict, "Can't do a strict left join."
        # No duplicates.
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
    else:  # inner join
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
        args.left_join, args.outer_join, args.unaligned_only,
        args.null_string)

    # Add the missing samples back to the matrix.
    if not args.dont_add_missing_samples:
        matrix_data = add_missing_samples(matrix_data, args.null_string)

    # Write out each of the matrices.
    for x in matrix_data:
        infile, outfile, matrix, header, samples = x
        if outfile == SKIP_OUTFILE:
            continue
        write_matrix(outfile, matrix)


if __name__ == '__main__':
    main()
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
