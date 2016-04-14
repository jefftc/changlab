#!/usr/bin/env python

# Functions:
# parse_indexes
#
# indexes_matrix
# select_cols_str
# select_cols_substr
# add_column
# copy_column
# 
# add_header_line
# fill_empty_headers
# remove_header_line
# reorder_headers_alphabetical
# upper_headers
# lower_headers
# hash_headers
# remove_duplicate_headers
# rename_duplicate_headers
# rename_header
# rename_header_i
# append_to_headers
# prepend_to_headers
# replace_header
# replace_header_re
#
# strip_all_annots
# upper_annots
# lower_annots
# set_value_if_empty
# copy_value_if_empty
# copy_value_if_empty_header
# copy_value_if_empty_same_header
# copy_value_if_empty_same_header_all
# replace_whole_annot
# replace_annots
# prepend_to_annots
# apply_re_to_annots
# merge_annots
# merge_annots_to_new_col
# split_annots
#
# _add_annots
# _subtract_annots
# _divide_annots
# _calc_two_annots
# 
# flip01_matrix
# all_same
# min_annots
# max_annots
# add_to
# multiply_by
# log_base
# neg_log_base
# add_two_annots
# subtract_two_annots
# divide_two_annots
# divide_many_annots
# average_same_header
# round_annots
#
# vcf_standardize
# vcf_remove_bad_coords
# vcf_remove_multicalls
# vcf_extract_format_values
# vcf_extract_info_values
# vcf_split_AD
# vcf_calc_vaf
#
# subtract_two_bed_lists
# subtract_value_from_bed_list


def parse_indexes(MATRIX, indexes_str, allow_duplicates=False,
                  check_range=True):
    # Takes 1-based indexes and returns a list of 0-based indexes.
    # 
    # Example inputs:
    # 5
    # 1,5,10
    # 1-99,215-300
    from genomicode import parselib

    max_index = len(MATRIX.headers)
    indexes_str = indexes_str.replace("END", str(max_index))

    I = []
    for s, e in parselib.parse_ranges(indexes_str):
        assert s >= 1
        if check_range:
            assert s <= len(MATRIX.headers), "Out of range: %d/%d" % (
                s, len(MATRIX.headers))
        s, e = s - 1, min(e, max_index)
        I.extend(range(s, e))

    if not allow_duplicates:
        # Remove duplicated indexes.  Need to preserve order.
        nodup = []
        for i in I:
            if i not in nodup:
                nodup.append(i)
        I = nodup

    return I


def indexes_matrix(MATRIX, indexes_list):
    # indexes is a list of strings indicating indexes.  Parse this and
    # return a submatrix consisting of just those indexes.
    from genomicode import AnnotationMatrix
    
    if not indexes_list:
        return MATRIX
    I = []
    for indexes in indexes_list:
        x = parse_indexes(MATRIX, indexes, allow_duplicates=True)
        I.extend(x)
    x = AnnotationMatrix.colslice(MATRIX, I)
    return x


def select_cols_str(MATRIX, cols_str):
    # cols_str is a list of the names of the headers to keep.
    if not cols_str:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    I = []
    for i, h in enumerate(MATRIX.headers):
        found = False
        for s in cols_str:
            if h == s:
                found = True
        if found:
            I.append(i)
    return AnnotationMatrix.colslice(MATRIX, I)


def select_cols_substr(MATRIX, cols_substr):
    # cols_substr is a list of the substrings of the headers to keep.
    if not cols_substr:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    I = []
    for i, h in enumerate(MATRIX.headers):
        found = False
        for s in cols_substr:
            if h.find(s) >= 0:
                found = True
        if found:
            I.append(i)
    return AnnotationMatrix.colslice(MATRIX, I)


def flip01_matrix(MATRIX, indexes):
    if not indexes:
        return MATRIX
    I = parse_indexes(MATRIX, indexes)

    MATRIX = MATRIX.copy()
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
        header_h = MATRIX.headers_h[i]
        annots = MATRIX.header2annots[header_h]
        for j in range(len(annots)):
            if annots[j].strip() == "0":
                annots[j] = "1"
            elif annots[j].strip() == "1":
                annots[j] = "0"
        MATRIX.header2annots[header_h] = annots
    return MATRIX


def fill_empty_headers(MATRIX, fill_headers):
    if not fill_headers:
        return MATRIX
    from genomicode import AnnotationMatrix

    headers = MATRIX.headers[:]
    for i in range(len(headers)):
        if headers[i].strip():
            continue
        j = 0
        while True:
            x = "H%03d" % j
            if x not in headers:
                break
            j += 1
        headers[i] = x
    return AnnotationMatrix.replace_headers(MATRIX, headers)


def reorder_headers_alphabetical(MATRIX, reorder_headers):
    if not reorder_headers:
        return MATRIX
    from genomicode import jmath
    from genomicode import AnnotationMatrix

    O = jmath.order(MATRIX.headers)
    headers = [MATRIX.headers[i] for i in O]
    headers_h = [MATRIX.headers_h[i] for i in O]
    M = AnnotationMatrix.AnnotationMatrix(
        headers, headers_h, MATRIX.header2annots)
    return M


def upper_headers(MATRIX, upper_headers):
    if not upper_headers:
        return MATRIX
    from genomicode import AnnotationMatrix

    # Convert to the upper case name.  Need to be careful because may
    # cause duplicates.
    headers = [x.upper() for x in MATRIX.headers]
    return AnnotationMatrix.replace_headers(MATRIX, headers)


def lower_headers(MATRIX, lower_headers):
    if not lower_headers:
        return MATRIX
    from genomicode import AnnotationMatrix

    # Convert to the lower case name.  Need to be careful because may
    # cause duplicates.
    headers = [x.lower() for x in MATRIX.headers]
    return AnnotationMatrix.replace_headers(MATRIX, headers)


def hash_headers(MATRIX, hash_headers):
    if not hash_headers:
        return MATRIX
    from genomicode import hashlib
    from genomicode import AnnotationMatrix

    # Hash each name.  Need to be careful because may cause
    # duplicates.
    headers = [hashlib.hash_var(x) for x in MATRIX.headers]
    return AnnotationMatrix.replace_headers(MATRIX, headers)


def add_header_line(filename, header_list, is_csv=False):
    # header_list is a list of a comma-separated list of headers.
    from genomicode import AnnotationMatrix
    from genomicode import filelib
    from genomicode import jmath

    delimiter = "\t"
    if is_csv:
        delimiter = ","

    X = [x for x in filelib.read_cols(filename, delimiter=delimiter)]
    # Check the dimensions of the matrix.
    assert X, "empty matrix"
    for i in range(len(X)):
        assert len(X[i]) == len(X[0])
    # Make each row an annotation.
    X = jmath.transpose(X)

    header_str = ",".join(header_list)
    x = header_str.split(",")
    assert len(x) >= len(X), "Matrix has %d columns, but %d headers given." % (
        len(X), len(x))
    # If there are more headers than columns, then fill the rest with
    # blanks.
    headers = x
    headers_h = AnnotationMatrix.uniquify_headers(headers)
    header2annots = {}
    for i, header_h in enumerate(headers_h):
        header_h = headers_h[i]
        annots = [""] * len(X[0])
        if i < len(X):
            annots = X[i]
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(headers, headers_h, header2annots)


def remove_header_line(filename, read_as_csv):
    from genomicode import AnnotationMatrix
    from genomicode import jmath

    MATRIX = AnnotationMatrix.read(filename, read_as_csv)
    matrix = []
    for header_h in MATRIX.headers_h:
        x = MATRIX.header2annots[header_h]
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)
    for x in matrix:
        print "\t".join(map(str, x))
    

def remove_duplicate_headers(MATRIX, remove_dups):
    if not remove_dups:
        return MATRIX
    from genomicode import AnnotationMatrix

    I = []
    seen = {}
    for i, h in enumerate(MATRIX.headers):
        if h in seen:
            continue
        seen[h] = 1
        I.append(i)
    return AnnotationMatrix.colslice(MATRIX, I)


def rename_duplicate_headers(MATRIX, rename_dups):
    if not rename_dups:
        return MATRIX
    from genomicode import AnnotationMatrix

    name2I = {}  # name -> list of indexes
    for i, name in enumerate(MATRIX.headers):
        if name not in name2I:
            name2I[name] = []
        name2I[name].append(i)

    nodup = MATRIX.headers[:]
    for (name, I) in name2I.iteritems():
        if len(I) < 2:
            continue
        for i in range(len(I)):
            nodup[I[i]] = "%s (%d)" % (name, i+1)

    x = AnnotationMatrix.replace_headers(MATRIX, nodup)
    return x
    

def rename_header(MATRIX, rename_list):
    # rename_list is list of strings in format of: <from>,<to>.
    if not rename_list:
        return MATRIX
    from genomicode import AnnotationMatrix

    rename_all = []  # list of (from_str, to_str)
    for rename_str in rename_list:
        x = rename_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        from_str, to_str = x
        rename_all.append((from_str, to_str))

    for from_str, to_str in rename_all:
        assert from_str in MATRIX.headers, "%s not a header" % from_str
        assert from_str in MATRIX.header2annots, "%s not a unique header" % \
               from_str

    convert = {}
    for from_str, to_str in rename_all:
        assert from_str not in convert, "dup: %s" % from_str
        convert[from_str] = to_str

    # Convert to the new names.
    headers = [convert.get(x, x) for x in MATRIX.headers]
    x = AnnotationMatrix.replace_headers(MATRIX, headers)
    return x
        

def rename_header_i(MATRIX, rename_list):
    # rename_list is list of strings in format of: <index>,<to>.
    if not rename_list:
        return MATRIX
    from genomicode import AnnotationMatrix

    rename_all = []  # list of (0-based index, to_str)
    for rename_str in rename_list:
        x = rename_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        index, to_str = x
        index = int(index)
        assert index >= 1 and index <= len(MATRIX.headers)
        index -= 1
        rename_all.append((index, to_str))

    # Convert to the new names.
    headers = MATRIX.headers[:]
    for index, to_str in rename_all:
        headers[index] = to_str
    x = AnnotationMatrix.replace_headers(MATRIX, headers)
    return x


def replace_header(MATRIX, replace_list):
    # replace_list is list of strings in format of: <from>,<to>.
    if not replace_list:
        return MATRIX
    from genomicode import AnnotationMatrix

    replace_all = []  # list of (from_str, to_str)
    for replace_str in replace_list:
        x = replace_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        from_str, to_str = x
        replace_all.append((from_str, to_str))

    # Convert to the new names.
    headers = MATRIX.headers[:]
    for from_str, to_str in replace_all:
        for i in range(len(headers)):
            x = headers[i]
            x = x.replace(from_str, to_str)
            headers[i] = x
    x = AnnotationMatrix.replace_headers(MATRIX, headers)
    return x
        

def replace_header_re(MATRIX, replace_list):
    # replace_list is list of strings in format of: <from re>,<to>.
    import re

    if not replace_list:
        return MATRIX
    from genomicode import AnnotationMatrix

    replace_all = []  # list of (from_re_str, to_str)
    for replace_str in replace_list:
        x = replace_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        from_str, to_str = x
        replace_all.append((from_str, to_str))

    # Convert to the new names.
    headers = MATRIX.headers[:]
    for from_str, to_str in replace_all:
        for i in range(len(headers)):
            x = headers[i]
            m = re.search(from_str, x)
            if m:
                x = x[:m.start(0)] + to_str + x[m.end(0):]
                #x = x.replace(m.group(0), to_str)
            headers[i] = x
    x = AnnotationMatrix.replace_headers(MATRIX, headers)
    return x


def append_to_headers(MATRIX, append_to_headers):
    # append_to_headers is list of strings in format of: <indexes>;<postfix>.
    if not append_to_headers:
        return MATRIX
    from genomicode import AnnotationMatrix

    append_all = []  # list of (list of 0-based indexes, postfix)
    for x in append_to_headers:
        x = x.split(";", 1)
        assert len(x) == 2
        indexes_str, prefix = x
        indexes = parse_indexes(MATRIX, indexes_str)
        for i in indexes:
            assert i >= 0 and i < len(MATRIX.headers)
        append_all.append((indexes, prefix))

    headers = MATRIX.headers[:]
    for indexes, postfix in append_all:
        for i in indexes:
            headers[i] = "%s%s" % (headers[i], prefix)
    return AnnotationMatrix.replace_headers(MATRIX, headers)


def prepend_to_headers(MATRIX, prepend_to_headers):
    # prepend_to_headers is list of strings in format of: <indexes>;<prefix>.
    if not prepend_to_headers:
        return MATRIX
    from genomicode import AnnotationMatrix

    prepend_all = []  # list of (list of 0-based indexes, prefix)
    for x in prepend_to_headers:
        x = x.split(";", 1)
        assert len(x) == 2
        indexes_str, prefix = x
        indexes = parse_indexes(MATRIX, indexes_str)
        for i in indexes:
            assert i >= 0 and i < len(MATRIX.headers)
        prepend_all.append((indexes, prefix))

    headers = MATRIX.headers[:]
    for indexes, prefix in prepend_all:
        for i in indexes:
            headers[i] = "%s%s" % (prefix, headers[i])
    x = AnnotationMatrix.replace_headers(MATRIX, headers)
    return x
    

def add_column(MATRIX, add_column):
    # add_column is list of strings in format of: <index>,<header>,<default>.
    if not add_column:
        return MATRIX
    from genomicode import AnnotationMatrix

    num_annots = None
    for annots in MATRIX.header2annots.itervalues():
        if num_annots is None:
            num_annots = len(annots)
        assert num_annots == len(annots)

    add_all = []  # list of (0-based index, header, default_value)
    last_index = -1
    for x in add_column:
        x = x.split(",", 2)
        assert len(x) == 3, "Format should be: <index>,<header>,<value>"
        index, header, default_value = x
        if index == "END":
            x = max(last_index+1, MATRIX.num_headers())
            index = x + 1
        index = int(index) - 1
        last_index = index
        add_all.append((index, header, default_value))
        
    # Since the hashed header names might change, keep track of the
    # indexes for each header.
    h_indexes = [("OLD", i) for i in range(len(MATRIX.headers))]
    for i, x in enumerate(add_all):
        index, header, default_value = x
        assert index >= 0 and index <= len(h_indexes)
        h_indexes.insert(index, ("NEW", i))

    headers = []
    for (which_one, i) in h_indexes:
        if which_one == "OLD":
            headers.append(MATRIX.headers[i])
        elif which_one == "NEW":
            index, header, default_value = add_all[i]
            headers.append(header)
        else:
            raise AssertionError
    headers_h = AnnotationMatrix.uniquify_headers(headers)

    header2annots = {}
    for i_new, (which_one, i_old) in enumerate(h_indexes):
        if which_one == "OLD":
            old_header_h = MATRIX.headers_h[i_old]
            new_header_h = headers_h[i_new]
            header2annots[new_header_h] = MATRIX.header2annots[old_header_h]
        elif which_one == "NEW":
            index, header, default_value = add_all[i_old]
            annots = [default_value] * num_annots
            new_header_h = headers_h[i_new]
            header2annots[new_header_h] = annots
        else:
            raise AssertionError

    return AnnotationMatrix.AnnotationMatrix(headers, headers_h, header2annots)


def copy_column(MATRIX, copy_column):
    # copy_column is a list of: <index>,<new_header>.
    if not copy_column:
        return MATRIX
    from genomicode import AnnotationMatrix

    jobs = []  # list of (0-based index, new_header)
    for x in copy_column:
        x = x.split(",", 1)
        assert len(x) == 2
        index, new_header = x
        index = int(index)
        assert index >= 1 and index <= len(MATRIX.headers)
        index -= 1  # convert to 0-based
        x = index, new_header
        jobs.append(x)

    headers = MATRIX.headers[:]
    all_annots = [MATRIX.header2annots[x] for x in MATRIX.headers_h]
    for x in jobs:
        index, new_header = x
        headers.append(new_header)
        all_annots.append(all_annots[index][:])

    headers_h = AnnotationMatrix.uniquify_headers(headers)
    assert len(headers_h) == len(all_annots)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(headers, headers_h, header2annots)


def set_value_if_empty(MATRIX, params):
    # list of strings in format of: <index 1-based>,<value>
    if not params:
        return MATRIX

    jobs = []    # list of (index 0-based, value)
    for x in params:
        x = x.split(",")
        assert len(x) == 2, "format should be: <index 1-based>,<value>"
        index, value = x
        index = int(index)
        # Should be 1-based indexes.
        assert index >= 1 and index <= len(MATRIX.headers)
        # Convert to 0-based indexes.
        index = index - 1
        jobs.append((index, value))

    MATRIX = MATRIX.copy()
    for x in jobs:
        index, value = x
        h = MATRIX.headers_h[index]
        
        # Change the annotations in place.
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            if not annots[i]:
                annots[i] = value
    return MATRIX


def copy_value_if_empty(MATRIX, copy_values):
    # copy_values is list of strings in format of: <dst>,<src 1>[,<src
    # 2>...].
    if not copy_values:
        return MATRIX

    copy_indexes = []   # list of (dst, src1 [, src 2...]).  0-based
    for copy_value in copy_values:
        x = copy_value.split(",")
        assert len(x) >= 2, "format should be: <dst>,<src 1>[, <src 2>...]"
        x = [int(x) for x in x]
        for i in range(len(x)):
            # Should be 1-based indexes.
            assert x[i] >= 1 and x[i] <= len(MATRIX.headers)
        # Convert to 0-based indexes.
        x = [x-1 for x in x]
        copy_indexes.append(tuple(x))

    MATRIX = MATRIX.copy()
    for indexes in copy_indexes:
        i_dst = indexes[0]
        header_dst = MATRIX.headers_h[i_dst]
        for i_src in indexes[1:]:
            header_src = MATRIX.headers_h[i_src]

            # Change the annotations in place.
            annots_dst = MATRIX.header2annots[header_dst]
            annots_src = MATRIX.header2annots[header_src]
            for i in range(len(annots_dst)):
                if not annots_dst[i].strip():
                    annots_dst[i] = annots_src[i]
    return MATRIX


def copy_value_if_empty_header(MATRIX, copy_values):
    # copy_values is list of strings in format of: <dst header>,<src
    # 1>[,<src 2>...].
    if not copy_values:
        return MATRIX

    copy_indexes = []   # list of (dst, src1 [, src 2...]).  0-based
    for copy_value in copy_values:
        headers = copy_value.split(",")
        assert len(headers) >= 2, \
               "format should be: <dst>,<src 1>[, <src 2>...]"
        indexes = []
        for header in headers:
            i = [i for i in range(len(MATRIX.headers))
                 if header == MATRIX.headers[i]]
            assert i, "Header not found: %s" % header
            assert len(i) == 1, "Header duplicated: %s" % header
            i = i[0]
            indexes.append(i)
        copy_indexes.append(tuple(indexes))

    MATRIX = MATRIX.copy()
    for indexes in copy_indexes:
        i_dst = indexes[0]
        header_dst = MATRIX.headers_h[i_dst]
        for i_src in indexes[1:]:
            header_src = MATRIX.headers_h[i_src]

            # Change the annotations in place.
            annots_dst = MATRIX.header2annots[header_dst]
            annots_src = MATRIX.header2annots[header_src]
            for i in range(len(annots_dst)):
                if not annots_dst[i].strip():
                    annots_dst[i] = annots_src[i]
    return MATRIX


def copy_value_if_empty_same_header(MATRIX, copy_values):
    # copy_values is list of header names.
    import itertools
    if not copy_values:
        return MATRIX

    copy_indexes = []   # list of list of indexes. 0-based
    for copy_value in copy_values:
        indexes = [i for i in range(len(MATRIX.headers))
                   if copy_value == MATRIX.headers[i]]
        assert indexes, "Header not found: %s" % copy_value
        assert len(indexes) > 1, "Header only found once: %s" % copy_value
        copy_indexes.append(indexes)

    MATRIX = MATRIX.copy()
    # Clean up the data.
    seen = {}
    all_indexes = []
    for I in copy_indexes:
        all_indexes.extend(I)
    all_indexes = {}.fromkeys(all_indexes)
    for i in all_indexes:
        header = MATRIX.headers_h[i]
        x = MATRIX.header2annots[header]
        x = [x.strip() for x in x]
        MATRIX.header2annots[header] = x
    
    for indexes in copy_indexes:
        all_headers = [MATRIX.headers_h[i] for i in indexes]
        all_annots = [MATRIX.header2annots[x] for x in all_headers]
        
        for i_dst, annots_dst in enumerate(all_annots):
            I = [i for (i, x) in enumerate(annots_dst) if not x]
            for k in I:
                for i_src, annots_src in enumerate(all_annots):
                    if i_src == i_dst:
                        continue
                    if annots_src[k]:
                        annots_dst[k] = annots_src[k]
                        break
        for header, annots in zip(all_headers, all_annots):
            MATRIX.header2annots[header] = annots
    return MATRIX


def copy_value_if_empty_same_header_all(MATRIX, copy_values):
    # copy_values is boolean
    if not copy_values:
        return MATRIX

    dup = []
    seen = {}
    for h in MATRIX.headers:
        if h in seen:
            dup.append(h)
        seen[h] = 1
    dup = {}.fromkeys(dup)
    return copy_value_if_empty_same_header(MATRIX, dup)


def strip_all_annots(MATRIX, strip):
    if not strip:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    header2annots = {}
    for header_h, annots in MATRIX.header2annots.iteritems():
        annots = [x.strip() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(
        MATRIX.headers, MATRIX.headers_h, header2annots)


def upper_annots(MATRIX, upper):
    if not upper:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    I = parse_indexes(MATRIX, upper)
    header2annots = MATRIX.header2annots.copy()
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
        header_h = MATRIX.headers_h[i]
        annots = MATRIX.header2annots[header_h]
        annots = [x.upper() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(
        MATRIX.headers, MATRIX.headers_h, header2annots)


def lower_annots(MATRIX, lower):
    if not lower:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    I = parse_indexes(MATRIX, lower)    
    header2annots = MATRIX.header2annots.copy()
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
        header_h = MATRIX.headers_h[i]
        annots = MATRIX.header2annots[header_h]
        annots = [x.lower() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(
        MATRIX.headers, MATRIX.headers_h, header2annots)


def replace_annot(MATRIX, replace_annot):
    # list of strings in format of: <indexes 1-based>;<src>;<dst>
    if not replace_annot:
        return MATRIX

    replace_all = []   # list of (indexes 0-based, src, dst)
    for replace in replace_annot:
        x = replace.split(";")
        assert len(x) == 3, "format should be: <indexes>;<src>;<dst>"
        indexes_str, src, dst = x
        indexes = parse_indexes(MATRIX, indexes_str)
        for index in indexes:
            replace_all.append((index, src, dst))

    MATRIX = MATRIX.copy()
    for x in replace_all:
        index, src, dst = x
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            # Change the annotations in place.
            annots[i] = annots[i].replace(src, dst)
    return MATRIX


def replace_whole_annot(MATRIX, replace_annot):
    # list of strings in format of: <indexes 1-based>;<src>;<dst>
    if not replace_annot:
        return MATRIX

    replace_all = []   # list of (indexes 0-based, src, dst)
    for replace in replace_annot:
        x = replace.split(";")
        assert len(x) == 3, "format should be: <indexes>;<src>;<dst>"
        indexes_str, src, dst = x
        indexes = parse_indexes(MATRIX, indexes_str)
        #index = int(index)
        ## Should be 1-based.
        #assert index >= 1 and index <= len(MATRIX.headers)
        ## Convert to 0-based.
        #index -= 1
        for index in indexes:
            replace_all.append((index, src, dst))

    MATRIX = MATRIX.copy()
    for x in replace_all:
        index, src, dst = x
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            # Change the annotations in place.
            if annots[i] == src:
                annots[i] = dst
    return MATRIX


def prepend_to_annots(MATRIX, prepend_annot):
    # list of strings in format of: <indexes 1-based>;<text to prepend>
    if not prepend_annot:
        return MATRIX

    prepend_all = []   # list of (index 0-based, src, dst)
    for prepend in prepend_annot:
        x = prepend.split(";")
        assert len(x) == 2, "format should be: <indexes>;<text>"
        indexes_str, text = x
        indexes = parse_indexes(MATRIX, indexes_str)
        for index in indexes:
            prepend_all.append((index, text))

    MATRIX = MATRIX.copy()
    for x in prepend_all:
        index, text = x
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            # Change the annotations in place.
            annots[i] = "%s%s" % (text, annots[i])
    return MATRIX


def apply_re_to_annots(MATRIX, apply_annots):
    # list of strings in format of: <indexes 1-based>;<regular expression>
    import re
    
    if not apply_annots:
        return MATRIX

    apply_all = []   # list of (index 0-based, regex)
    for apply_ in apply_annots:
        x = apply_.split(";")
        assert len(x) == 2, "format should be: <indexes>;<regex>"
        indexes_str, regex = x
        indexes = parse_indexes(MATRIX, indexes_str)
        for index in indexes:
            apply_all.append((index, regex))

    MATRIX = MATRIX.copy()
    for x in apply_all:
        index, regex = x
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            # Change the annotations in place.
            m = re.search(regex, annots[i])
            if m:
                annots[i] = m.group(1)
    return MATRIX


def merge_annots(MATRIX, merge_annots):
    # list of strings in format of:
    # <src indexes 1-based>;<dst index 1-based>;<char>
    if not merge_annots:
        return MATRIX

    merge_all = []   # list of (src indexes 0-based, dst index 0-based, char)
    for merge in merge_annots:
        x = merge.split(";")
        assert len(x) == 3, \
               "format should be: <src indexes>;<dst index>;<char>"
        src_indexes_str, dst_indexes_str, merge_char = x
        src_indexes = parse_indexes(MATRIX, src_indexes_str)
        dst_indexes = parse_indexes(MATRIX, dst_indexes_str)
        assert len(dst_indexes) == 1
        dst_index = dst_indexes[0]
        merge_all.append((src_indexes, dst_index, merge_char))

    MATRIX = MATRIX.copy()
    for x in merge_all:
        src_indexes, dst_index, merge_char = x

        src_annots = []
        for i in src_indexes:
            h = MATRIX.headers_h[i]
            x = MATRIX.header2annots[h]
            src_annots.append(x)
        # Change MATRIX place.
        h = MATRIX.headers_h[dst_index]
        dst_annots = MATRIX.header2annots[h]
        for i in range(len(dst_annots)):
            x = [x[i] for x in src_annots]
            merged = merge_char.join(x)
            dst_annots[i] = merged
    return MATRIX


def merge_annots_to_new_col(MATRIX, merge_annots):
    # list of strings in format of:
    # <src indexes 1-based>;<dst col name>;<char>
    if not merge_annots:
        return MATRIX
    from genomicode import AnnotationMatrix

    jobs = []   # list of (src indexes 0-based, dst_name, char)
    for fmt in merge_annots:
        x = fmt.split(";")
        assert len(x) == 3, \
               "format should be: <src indexes>;<dst name>;<char>.  Got %s" % \
               fmt
        src_indexes_str, dst_name, merge_char = x
        src_indexes = parse_indexes(MATRIX, src_indexes_str)
        jobs.append((src_indexes, dst_name, merge_char))

    headers = MATRIX.headers[:]
    all_annots = [MATRIX.header2annots[x] for x in MATRIX.headers_h]
    
    for x in jobs:
        src_indexes, dst_name, merge_char = x

        src_annots = []
        for i in src_indexes:
            h = MATRIX.headers_h[i]
            x = MATRIX.header2annots[h]
            src_annots.append(x)
        dst_annots = [""] * MATRIX.num_annots()
        for i in range(len(dst_annots)):
            x = [x[i] for x in src_annots]
            x = merge_char.join(x)
            dst_annots[i] = x
        headers.append(dst_name)
        all_annots.append(dst_annots)

    return AnnotationMatrix.create_from_annotations(headers, all_annots)


def split_annots(MATRIX, split_annots):
    # list of strings in format of:
    # <src index>;<dst indexes>;<split char>
    if not split_annots:
        return MATRIX

    jobs = []   # list of (src index 0-based, dst indexes 0-based, char)
    for x in split_annots:
        x = x.split(";")
        assert len(x) == 3, \
               "format should be: <src index>;<dst indexes>;<char>"
        src_index_str, dst_indexes_str, split_char = x
        src_indexes = parse_indexes(MATRIX, src_index_str)
        dst_indexes = parse_indexes(MATRIX, dst_indexes_str)
        assert len(src_indexes) == 1
        src_index = src_indexes[0]
        jobs.append((src_index, dst_indexes, split_char))

    MATRIX = MATRIX.copy()
    for x in jobs:
        src_index, dst_indexes, split_char = x

        h = MATRIX.headers_h[src_index]
        src_annots = MATRIX.header2annots[h]
        split_annots = [x.split(split_char) for x in src_annots]
        for i in range(len(split_annots)):
            assert len(split_annots[i]) == len(dst_indexes), \
                   "split/dst_indexes mismatch: %d %s %s" % (
                i, split_annots[i], len(dst_indexes))
        for i in range(len(dst_indexes)):
            h = MATRIX.headers_h[dst_indexes[i]]
            dst_annots = MATRIX.header2annots[h]
            assert len(split_annots) == len(dst_annots)
            for j in range(len(split_annots)):
                # change in place
                dst_annots[j] = split_annots[j][i]
    return MATRIX


def split_chr_start_end(MATRIX, arg):
    # list of strings in format of: <header>
    if not arg:
        return MATRIX
    from genomicode import AnnotationMatrix

    jobs = []   # list of (index 0-based,)
    for x in arg:
        h = MATRIX.normalize_header(x, index_base1=True)
        assert h is not None, "Unknown header: %s" % x
        i = MATRIX.headers_h.index(h)
        jobs.append((i,))


    headers = MATRIX.headers[:]
    all_annots = [MATRIX.header2annots[x] for x in MATRIX.headers_h]
    for x in jobs:
        index, = x

        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        all_chrom = [""] * len(annots)
        all_start = [""] * len(annots)
        all_end = [""] * len(annots)
        for i, x in enumerate(annots):
            # chr1:320117-320142
            x = x.strip()
            if not x:
                continue
            x = x.split(":")
            assert len(x) == 2, "Bad format: %s" % annots[i]
            chrom, pos = x
            x = pos.split("-")
            assert len(x) == 2, "Bad format: %s" % annots[i]
            start, end = x
            start, end = int(start), int(end)
            assert end >= start
            all_chrom[i] = chrom
            all_start[i] = start
            all_end[i] = end

        h = MATRIX.headers[index]
        x1 = "%s chr" % h
        x2 = "%s start" % h
        x3 = "%s end" % h
        headers.extend([x1, x2, x3])
        all_annots.extend([all_chrom, all_start, all_end])

    return AnnotationMatrix.create_from_annotations(headers, all_annots)


def _add_annots(a1, a2):
    return a1 + a2

def _subtract_annots(a1, a2):
    return a1 - a2

def _divide_annots(a1, a2):
    num, den = a1, a2
    if abs(den) < 1E-50:
        return ""
    return num / den


def _calc_two_annots(MATRIX, calc_annots, calc_fn):
    # calc_annots is a list of <annot 1>,<annot 2>,<dest>.  Each are
    # 1-based indexes.  Returns a Matrix with the calculation applied.
    if not calc_annots:
        return MATRIX

    to_calc = []  # list of (i1, i2, i_dest); 0-based
    for ca in calc_annots:
        x = ca.split(",")
        assert len(x) == 3, "format should be: <annot1>,<annot2>,<dest>"
        i_1, i_2, i_dest = x
        i_1, i_2, i_dest = int(i_1), int(i_2), int(i_dest)
        # Convert to 0-based index.
        i_1, i_2, i_dest = i_1-1, i_2-1, i_dest-1
        assert i_1 >= 0 and i_1 < len(MATRIX.headers)
        assert i_2 >= 0 and i_2 < len(MATRIX.headers)
        assert i_dest >= 0 and i_dest < len(MATRIX.headers)
        x = i_1, i_2, i_dest
        to_calc.append(x)

    MATRIX = MATRIX.copy()
    for (i_1, i_2, i_dest) in to_calc:
        h_1 = MATRIX.headers_h[i_1]
        h_2 = MATRIX.headers_h[i_2]
        h_dest = MATRIX.headers_h[i_dest]
        annots_1 = MATRIX.header2annots[h_1]
        annots_2 = MATRIX.header2annots[h_2]
        assert len(annots_1) == len(annots_2)
        annots_dest = [""] * len(annots_1)
        for i in range(len(annots_1)):
            a1 = annots_1[i]
            a2 = annots_2[i]
            if not a1.strip() or not a2.strip():
                continue
            a1 = float(a1)
            a2 = float(a2)
            annots_dest[i] = calc_fn(a1, a2)
        MATRIX.header2annots[h_dest] = annots_dest
    
    return MATRIX


def all_same(MATRIX, all_same):
    # format: <indexes 1-based>;<dest index>
    if not all_same:
        return MATRIX

    x = all_same.split(";")
    assert len(x) == 2, "format should be: <indexes>;<index dest>"
    indexes_str, dst_i = x
    indexes = parse_indexes(MATRIX, indexes_str)
    dst_i = int(dst_i)
    assert dst_i >= 1 and dst_i <= MATRIX.num_headers()
    dst_i -= 1

    MATRIX = MATRIX.copy()
    annot_matrix = []  # indexes x annot matrix
    for i in indexes:
        h = MATRIX.headers_h[i]
        x = MATRIX.header2annots[h]
        annot_matrix.append(x)

    same_annot = [1] * MATRIX.num_annots()
    for i in range(MATRIX.num_annots()):
        # See if all annot_matrix[i] is same.
        same = True
        for j in range(1, len(annot_matrix)):
            if annot_matrix[j][i] != annot_matrix[0][i]:
                same = False
        if not same:
            same_annot[i] = 0

    h = MATRIX.headers_h[dst_i]
    MATRIX.header2annots[h] = same_annot
    return MATRIX


def min_annots(MATRIX, min_annots):
    # format: <indexes 1-based>;<dest index>
    from genomicode import jmath
    
    if not min_annots:
        return MATRIX

    x = min_annots.split(";")
    assert len(x) == 2, "format should be: <indexes>;<index dest>"
    indexes_str, dst_i = x
    indexes = parse_indexes(MATRIX, indexes_str)
    dst_i = int(dst_i)
    assert dst_i >= 1 and dst_i <= MATRIX.num_headers()
    dst_i -= 1

    MATRIX = MATRIX.copy()
    annot_matrix = []  # indexes x annot matrix
    for i in indexes:
        h = MATRIX.headers_h[i]
        x = MATRIX.header2annots[h]
        x = map(float, x)
        annot_matrix.append(x)
    mins = jmath.min(annot_matrix, byrow=False)
    assert len(mins) == MATRIX.num_annots()

    h = MATRIX.headers_h[dst_i]
    MATRIX.header2annots[h] = mins
    return MATRIX


def max_annots(MATRIX, max_annots):
    # format: <indexes 1-based>;<dest index>
    from genomicode import jmath
    
    if not max_annots:
        return MATRIX

    x = max_annots.split(";")
    assert len(x) == 2, "format should be: <indexes>;<index dest>"
    indexes_str, dst_i = x
    indexes = parse_indexes(MATRIX, indexes_str)
    dst_i = int(dst_i)
    assert dst_i >= 1 and dst_i <= MATRIX.num_headers()
    dst_i -= 1

    MATRIX = MATRIX.copy()
    annot_matrix = []  # indexes x annot matrix
    for i in indexes:
        h = MATRIX.headers_h[i]
        x = MATRIX.header2annots[h]
        x = map(float, x)
        annot_matrix.append(x)
    maxes = jmath.max(annot_matrix, byrow=False)
    assert len(maxes) == MATRIX.num_annots()

    h = MATRIX.headers_h[dst_i]
    MATRIX.header2annots[h] = maxes
    return MATRIX


def add_to(MATRIX, add_to):
    # format: list of <header or index 1-based>,<number>
    if not add_to:
        return MATRIX
    
    jobs = []  # list of (0-based index, number)
    for x in add_to:
        x = x.split(",")
        assert len(x) == 2, "format should be: <index>,<number>"
        header, number = x
        index = MATRIX.normalize_header_i(header, index_base1=True)
        assert index is not None, "Unknown header or index: %s" % header
        number = _int_or_float(number)
        x = index, number
        jobs.append(x)
    
    MATRIX = MATRIX.copy()
    for x in jobs:
        index, number = x
        assert index < len(MATRIX.headers_h)
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            x = _int_or_float(annots[i])
            x = x + number
            annots[i] = str(x)
    return MATRIX


def multiply_by(MATRIX, multiply_by):
    # format: list of <index 1-based>,<number>
    if not multiply_by:
        return MATRIX
    
    jobs = []  # list of (0-based index, number)
    for x in multiply_by:
        x = x.split(",")
        assert len(x) == 2, "format should be: <index>,<number>"
        index, number = x
        index = int(index)
        number = float(number)
        x = index-1, number
        jobs.append(x)
    
    MATRIX = MATRIX.copy()
    for x in jobs:
        index, number = x
        assert index < len(MATRIX.headers_h)
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            x = float(annots[i])
            x = x * number
            annots[i] = str(x)
    return MATRIX


def log_base(MATRIX, log_base):
    # format: list of <index 1-based>,<base>
    if not log_base:
        return MATRIX
    import math
    
    jobs = []  # list of (0-based index, base)
    for x in log_base:
        x = x.split(",")
        assert len(x) == 2, "format should be: <index>,<base>"
        index, base = x
        index = int(index)
        base = float(base)
        x = index-1, base
        jobs.append(x)

    MIN = 1E-100
    MATRIX = MATRIX.copy()
    for x in jobs:
        index, base = x
        assert index < len(MATRIX.headers_h)
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            x = float(annots[i])
            x = max(x, MIN)
            x = math.log(x, base)
            annots[i] = str(x)
    return MATRIX


def neg_log_base(MATRIX, log_base):
    # format: list of <index 1-based>,<base>
    if not log_base:
        return MATRIX
    import math
    
    jobs = []  # list of (0-based index, base)
    for x in log_base:
        x = x.split(",")
        assert len(x) == 2, "format should be: <index>,<base>"
        index, base = x
        index = int(index)
        base = float(base)
        x = index-1, base
        jobs.append(x)

    MIN = 1E-100
    MATRIX = MATRIX.copy()
    for x in jobs:
        index, base = x
        assert index < len(MATRIX.headers_h)
        h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[h]
        for i in range(len(annots)):
            x = float(annots[i])
            x = max(x, MIN)
            x = math.log(x, base)
            x = x * -1
            annots[i] = str(x)
    return MATRIX


def add_two_annots(MATRIX, add_annots):
    # Format: list of <annot 1>,<annot 2>,<dest>.  Each are 1-based
    # indexes.  dest is annot1 - annot2.
    return _calc_two_annots(MATRIX, add_annots, _add_annots)
    

def subtract_two_annots(MATRIX, subtract_annots):
    # Format: list of <annot 1>,<annot 2>,<dest>.  Each are 1-based
    # indexes.  dest is annot1 - annot2.
    return _calc_two_annots(MATRIX, subtract_annots, _subtract_annots)


def divide_two_annots(MATRIX, divide_annots):
    # Format: list of <numerator>,<denominator>,<dest>.  Each are 1-based
    # indexes.
    return _calc_two_annots(MATRIX, divide_annots, _divide_annots)


def divide_many_annots(MATRIX, divide_annots):
    # Format: list of <numerator indexes>;<denominator index>.  Each
    # are 1-based indexes.
    if not divide_annots:
        return MATRIX
    divide_all = []  # list of (list of 0-based indexes, 0-based index)
    for x in divide_annots:
        x = x.split(";")
        assert len(x) == 2
        x1, x2 = x
        indexes1 = parse_indexes(MATRIX, x1)
        indexes2 = parse_indexes(MATRIX, x2)
        assert len(indexes2) == 1
        for i in indexes1 + indexes2:
            assert i >= 0 and i < len(MATRIX.headers)
        divide_all.append((indexes1, indexes2[0]))

    MATRIX = MATRIX.copy()
    for x in divide_all:
        num_indexes, den_index = x
        for i in range(MATRIX.num_annots()):
            header_h = MATRIX.headers_h[den_index]
            annots = MATRIX.header2annots[header_h]
            den = float(annots[i])

            for index in num_indexes:
                header_h = MATRIX.headers_h[index]
                annots = MATRIX.header2annots[header_h]
                num = float(annots[i])
                annots[i] = num / den
    
    return MATRIX


def average_same_header(MATRIX, average):
    if not average:
        return MATRIX
    from genomicode import jmath
    from genomicode import AnnotationMatrix

    # Make a list of all the duplicate headers.
    header2I = {}  # header -> list of indexes
    for i, header in enumerate(MATRIX.headers):
        if header not in header2I:
            header2I[header] = []
        header2I[header].append(i)

    # Now make the new matrix.
    headers = []
    header2annots = {}
    for header in MATRIX.headers:
        if header in header2annots:
            continue
        I = header2I[header]
        MATRIX_I = []
        for i in I:
            h = MATRIX.headers_h[i]
            x = MATRIX.header2annots[h]
            MATRIX_I.append(x)
        if len(MATRIX_I) == 1:
            x = MATRIX_I[0]
        else:
            for i in range(len(MATRIX_I)):
                MATRIX_I[i] = [float(x) for x in MATRIX_I[i]]
            x = jmath.mean(MATRIX_I, byrow=0)
        headers.append(header)
        header2annots[header] = x
    return AnnotationMatrix.AnnotationMatrix(headers, headers, header2annots)


def round_annots(MATRIX, round_annots):
    # Format: list of <index>.  1-based indexes.
    if not round_annots:
        return MATRIX
    
    indexes = []  # list of 0-based indexes
    for x in round_annots:
        I = parse_indexes(MATRIX, x)
        for i in I:
            assert i >= 0 and i < len(MATRIX.headers)
        indexes.extend(I)
    indexes = sorted({}.fromkeys(indexes))

    MATRIX = MATRIX.copy()
    for index in indexes:
        header_h = MATRIX.headers_h[index]
        x = MATRIX.header2annots[header_h]
        x = [int(round(float(x))) for x in x]
        MATRIX.header2annots[header_h] = x
    
    return MATRIX


def convert_percent_to_decimal(MATRIX, convert):
    # Format: list of <index>.  1-based indexes.
    if not convert:
        return MATRIX
    
    indexes = []  # list of 0-based indexes
    for x in convert:
        I = parse_indexes(MATRIX, x)
        for i in I:
            assert i >= 0 and i < len(MATRIX.headers)
        indexes.extend(I)
    indexes = sorted({}.fromkeys(indexes))

    MATRIX = MATRIX.copy()
    for index in indexes:
        header_h = MATRIX.headers_h[index]
        annots = MATRIX.header2annots[header_h]
        for i in range(len(annots)):
            x = annots[i]
            x = x.strip()
            if not x:
                continue
            if x.endswith("%"):
                x = x[:-1]
            x = float(x) / 100
            annots[i] = x
        MATRIX.header2annots[header_h] = annots
    return MATRIX


## def _header_or_index(MATRIX, header):
##     # header may be either a header or a 1-based index.  Return the
##     # hashed header.
##     if not header:
##         return None
##     if header in MATRIX.headers:
##         i = MATRIX.headers.index(header)
##         return MATRIX.headers_h[i]
##     if header in MATRIX.headers_h:
##         return header
##     header_i = None
##     try:
##         header_i = int(header)
##     except ValueError, x:
##         pass
##     if header_i is not None:
##         assert header_i >= 1 and header_i <= len(MATRIX.headers)
##         i = header_i - 1
##         return MATRIX.headers_h[i]
##     raise AssertionError, "Unknown header: %s" % header


def vcf_standardize(MATRIX, vcf_standardize):
    if not vcf_standardize:
        return MATRIX
    from genomicode import AnnotationMatrix
    from genomicode import vcflib

    # Format: <info_header>,<format_header>[,<genotype_header>]

    x = vcf_standardize.split(",")
    assert len(x) >= 2, \
           "Format: <info_header>,<format_header>,<genotype_header>"
    info_header, format_header = x[:2]
    genotype_headers = x[2:]

    info_header_n = MATRIX.normalize_header(info_header)
    format_header_n = MATRIX.normalize_header(format_header)
    assert info_header_n, "Missing header: %s" % info_header
    assert format_header_n, "Missing header: %s" % format_header

    if not genotype_headers:
        # Find the genotype headers at the end of the file.
        i1 = MATRIX.headers_h.index(info_header_n)
        i2 = MATRIX.headers_h.index(format_header_n)
        i_start = max(i1, i2) + 1
        assert i_start < len(MATRIX.headers), "No columns at end of file."
        for i in range(len(MATRIX.headers)-1, i_start-1, -1):
            # See if every row is either blank or contains some colons.
            h = MATRIX.headers_h[i]
            annots = MATRIX.header2annots[h]
            x = [x for x in annots if not x.strip() or x.find(":") >= 0]
            if len(x) != len(annots):
                break
        i_geno = i+1
        genotype_headers = MATRIX.headers[i_geno:]
    assert genotype_headers, "No genotype headers found."

    # Create a VCF object.
    samples = genotype_headers
    # Parse the info line.
    x = MATRIX.header2annots[info_header_n]
    more_info = [vcflib._parse_info_dict(x) for x in x]
    # Parse the genotype data.
    format_strings = MATRIX.header2annots[format_header_n]
    genotypes = {}
    for sample in genotype_headers:
        genotype_strings = MATRIX[sample]
        geno_dicts = [
            vcflib._parse_genotype_dict(fs, gs)
            for (fs, gs) in zip(format_strings, genotype_strings)]
        genotypes[sample] = geno_dicts
    vcf = vcflib.VCFFile(MATRIX, samples, more_info, genotypes)
    

    CHROM = "chrom"
    START = "start"
    END = "end"
    GENE = "gene"
    GENE_ID = "entrez_gene_id"
    FUNC = "func"
    EXONICFUNC = "exonicfunc"
    AACHANGE = "aachange"
    NUM_REF = "num_ref"
    NUM_ALT = "num_alt"
    TOTAL = "total_reads"
    VAF = "vaf"
    CALL = "call"

    # If I can't find these, then just fill with blank spaces.
    # This can happen if the file is not annotated.
    IGNORE_IF_MISSING = [GENE, GENE_ID, FUNC, EXONICFUNC, AACHANGE]

    # List of tuples:
    # - header name
    # - list of possible original headers
    COMMON_COLUMNS = [
        (CHROM, ["chrom", "contig", "Chr", "CHROM", "#CHROM"]),
        (START, ["start", "position", "Start", "POS", "pos"]),
        (END, ["end", "End"]),
        ("ref_allele", ["ref_allele", "Ref", "REF"]),
        ("alt_allele", ["alt_allele", "Alt", "ALT"]),
        (GENE, ["gene", "Gene", "Gene.refGene"]),
        (GENE_ID, ["entrez_gene_id", "Entrez_Gene_Id"]),
        (FUNC, ["func", "Func", "Func.refGene"]),
        (EXONICFUNC, ["exonicfunc", "ExonicFunc", "ExonicFunc.refGene"]),
        (AACHANGE, ["aachange", "AAChange", "AAChange.refGene"]),
        ]
    # Sample specific columns.
    SPECIFIC_COLUMNS = [
        (NUM_REF, ["num_ref", "t_ref_count"], "num_ref"),
        (NUM_ALT, ["num_alt", "t_alt_count"], "num_alt"),
        (TOTAL, ["total_reads"], "total_reads"),
        (VAF, ["vaf"], "vaf"),
        (CALL, [], "call"),
        ]

    headers = []
    header2annots = {}  # should contain no duplicates
    missing = []
 
    # Set the common columns.
    for (dst_header, src_headers) in COMMON_COLUMNS:
        header_i = None
        for h in src_headers:
            if h in MATRIX.headers:
                header_i = MATRIX.headers.index(h)
                break
        headers.append(dst_header)
        if header_i is None:
            missing.append(dst_header)
            continue
        assert dst_header not in header2annots
        h = MATRIX.headers_h[header_i]
        annots = MATRIX.header2annots[h]
        header2annots[dst_header] = annots

    # Set the sample-specific columns.
    for sample in genotype_headers:
        info_list = [
            vcflib.parse_info(vcf, sample, i)
            for i in range(MATRIX.num_annots())]
        for (dst_header, src_headers, info_member) in SPECIFIC_COLUMNS:
            if len(genotype_headers) > 1:
                dst_header = "%s %s" % (sample, dst_header)
            
            # If there is only one sample, look for the src_headers.
            if len(genotype_headers) == 1:
                header_i = None
                for h in src_headers:
                    if h in MATRIX.headers:
                        header_i = MATRIX.headers.index(h)
                        break
                if header_i:
                    headers.append(dst_header)
                    assert dst_header not in header2annots
                    h = MATRIX.headers_h[header_i]
                    annots = MATRIX.header2annots[h]
                    header2annots[dst_header] = annots
                    continue
            # Pull the information out the the info_list.
            headers.append(dst_header)
            assert dst_header not in header2annots
            x = [getattr(x, info_member) for x in info_list]
            x = [vcflib._fmt_vcf_value(x) for x in x]
            header2annots[dst_header] = x

    # If I can't find "end", and I could find the "start", then make
    # it the same as start.
    assert START not in missing
    if END in missing:
        header2annots[END] = header2annots[START][:]
        missing.pop(missing.index(END))

    # Ignore missing headers.
    for header in IGNORE_IF_MISSING:
        if header not in missing:
            continue
        missing.pop(missing.index(header))
        annots = [""] * MATRIX.num_annots()
        header2annots[header] = annots

    # Make sure nothing is missing.
    assert not missing, "Not found: %s" % ", ".join(map(str, missing))

    all_annots = [header2annots.get(x) for x in headers]

    # Clean up all annots.
    for i in range(len(all_annots)):
        for j in range(len(all_annots[i])):
            if all_annots[i][j] is None:
                all_annots[i][j] = ""
            all_annots[i][j] = str(all_annots[i][j]).strip()

    return AnnotationMatrix.create_from_annotations(headers, all_annots)


def vcf_remove_bad_coords(MATRIX, vcf_remove_bad_coords):
    if not vcf_remove_bad_coords:
        return MATRIX
    from genomicode import AnnotationMatrix

    # Column names must be standardized.
    START = "start"
    END = "end"
    assert START in MATRIX.headers_h, "VCF must have standardized names"
    assert END in MATRIX.headers_h, "VCF must have standardized names"

    start_annots = MATRIX.header2annots[START]
    end_annots = MATRIX.header2annots[END]
    start_annots = [x.lower() for x in start_annots]
    end_annots = [x.lower() for x in end_annots]
    
    bad_indexes = {}
    for i in range(len(start_annots)):
        s, e = start_annots[i], end_annots[i]
        if s.find("e") >= 0:
            bad_indexes[i] = 1
        elif e.find("e") >= 0:
            bad_indexes[i] = 1
    if not bad_indexes:
        return MATRIX

    MATRIX = MATRIX.copy()
    for h, annots in MATRIX.header2annots.iteritems():
        annots = [x for (i, x) in enumerate(annots) if i not in bad_indexes]
        MATRIX.header2annots[h] = annots
    return MATRIX


def vcf_remove_multicalls(MATRIX, vcf_remove_multicalls):
    if not vcf_remove_multicalls:
        return MATRIX
    from genomicode import AnnotationMatrix

    # Column names must be standardized.
    CHROM = "chrom"
    START = "start"
    END = "end"
    TOTAL = "total_reads"
    assert CHROM in MATRIX.headers_h, "VCF must have standardized names"
    assert START in MATRIX.headers_h, "VCF must have standardized names"
    assert END in MATRIX.headers_h, "VCF must have standardized names"

    chrom_annots = MATRIX.header2annots[CHROM]
    start_annots = MATRIX.header2annots[START]
    end_annots = MATRIX.header2annots[END]
    total_annots = MATRIX.header2annots[TOTAL]

    start_annots = [int(x) for x in start_annots]
    end_annots = [int(x) for x in end_annots]

    # Find the duplicates.
    loc2indexes = {}  # (chrom, start, end) -> list of indexes
    for i in range(len(chrom_annots)):
        chrom, start, end = chrom_annots[i], start_annots[i], end_annots[i]
        x = chrom, start, end
        if x not in loc2indexes:
            loc2indexes[x] = []
        loc2indexes[x].append(i)

    # Find rows to discard.
    bad_indexes = {}
    for (loc, indexes) in loc2indexes.iteritems():
        if len(indexes) < 2:
            continue
        # If there are duplicates, choose the best one.
        most_i = most_reads = None
        for i in indexes:
            if most_reads is None or total_annots[i] > most_reads:
                most_reads = total_annots[i]
                most_i = i
        assert most_i is not None
        for i in indexes:
            if i != most_i:
                bad_indexes[i] = 1

    MATRIX = MATRIX.copy()
    # Edit MATRIX in place.
    for h, annots in MATRIX.header2annots.iteritems():
        annots = [x for (i, x) in enumerate(annots) if i not in bad_indexes]
        MATRIX.header2annots[h] = annots
    return MATRIX
    

def vcf_extract_format_values(MATRIX, vcf_format):
    # Format: <format header>,<values header>,<value>[,value].
    if not vcf_format:
        return MATRIX
    from genomicode import AnnotationMatrix

    x = vcf_format.split(",")
    x = [x.strip() for x in x]
    assert len(x) >= 3, "Format: <header>,<header>,<value>[,<value>...]"

    f_header = x[0]
    v_header = x[1]
    value_headers = x[2:]

    assert f_header in MATRIX.headers, "Missing header: %s" % f_header
    assert v_header in MATRIX.headers, "Missing header: %s" % v_header
    # Assume no duplicates.  Just use the first one.

    h_f = MATRIX.headers_h[MATRIX.headers.index(f_header)]
    h_v = MATRIX.headers_h[MATRIX.headers.index(v_header)]
    annots_f = MATRIX.header2annots[h_f]  # list of strings
    annots_v = MATRIX.header2annots[h_v]  # list of strings
    assert len(annots_f) == len(annots_v)

    # Parse out the annotations into a matrix.
    annots_f = [x.split(":") for x in annots_f]
    annots_v = [x.split(":") for x in annots_v]

    headers = MATRIX.headers[:]
    all_annots = [MATRIX.header2annots[x] for x in MATRIX.headers_h]
    for value_header in value_headers:
        values = [""] * len(annots_f)
        for i in range(len(annots_f)):
            fmt = annots_f[i]
            vals = annots_v[i]
            # 1.  Sometimes len(vals) < len(fmt).
            #     GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
            #     ./.:.:1
            # 2.  Sometimes the value_header is missing in one line.
            #     GT:GQ:PL  (but AD in every other line)
            #assert len(fmt) == len(vals)
            if value_header not in fmt:
                continue
            #assert value_header in fmt, \
            #       "Missing value for: %s %s" % (value_header, annots_f[i])
            j = fmt.index(value_header)
            if j < len(vals):
                values[i] = vals[j]
        headers.append(value_header)
        all_annots.append(values)
                       
    headers_h = AnnotationMatrix.uniquify_headers(headers)
    assert len(headers_h) == len(all_annots)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    return AnnotationMatrix.AnnotationMatrix(headers, headers_h, header2annots)


def vcf_extract_info_values(MATRIX, vcf_info):
    # Format: <format header>,<value>[,value].
    if not vcf_info:
        return MATRIX
    from genomicode import AnnotationMatrix

    x = vcf_info.split(",")
    x = [x.strip() for x in x]
    assert len(x) >= 2, "Format: <header>,<value>[,<value>...]"

    i_header = x[0]
    value_names = x[1:]

    assert i_header in MATRIX.headers, "Missing header: %s" % i_header
    # Assume no duplicates.  Just use the first one.

    # Parse out the annotations.
    h_i = MATRIX.headers_h[MATRIX.headers.index(i_header)]
    # list of strings.  <name>=<value>[;<name>=<value>]
    annots_str = MATRIX.header2annots[h_i]  # list of strings
    # list of dicts.  <name>=<value>
    annots_dict = []  # list of dicts
    for x in annots_str:
        x = x.split(";")
        d = {}
        for x in x:
            x = x.split("=")
            assert len(x) == 2
            key, value = x
            d[key] = value
        annots_dict.append(d)
        
    headers = MATRIX.headers[:]
    all_annots = [MATRIX.header2annots[x] for x in MATRIX.headers_h]
    for value_name in value_names:
        values = [d.get(value_name, "") for d in annots_dict]
        headers.append(value_name)
        all_annots.append(values)
        
    return AnnotationMatrix.create_from_annotations(headers, all_annots)


def vcf_split_AD(MATRIX, split_annots):
    # list of strings in format of:
    # <src index>;<dst indexes>
    if not split_annots:
        return MATRIX

    jobs = []   # list of (src index 0-based, dst indexes 0-based, char)
    for x in split_annots:
        x = x.split(";")
        assert len(x) == 2, \
               "format should be: <src index>;<dst indexes>"
        src_index_str, dst_indexes_str = x
        src_indexes = parse_indexes(MATRIX, src_index_str)
        dst_indexes = parse_indexes(MATRIX, dst_indexes_str)
        assert len(src_indexes) == 1
        src_index = src_indexes[0]
        split_char = ","
        jobs.append((src_index, dst_indexes, split_char))

    MATRIX = MATRIX.copy()
    for x in jobs:
        src_index, dst_indexes, split_char = x

        h = MATRIX.headers_h[src_index]
        src_annots = MATRIX.header2annots[h]
        split_annots = [x.split(split_char) for x in src_annots]
        for i in range(len(split_annots)):
            if len(split_annots[i]) == len(dst_indexes):
                continue
            # If there are only 2 dst_indexes, they should refer REF
            # and ALT alleles.  In this case, just add everything up
            # into the ALT allele.
            if len(dst_indexes) == 2:
                x0 = split_annots[i][0]
                x1 = sum(map(int, split_annots[i][1:]))
                split_annots[i] = [x0, x1]
            assert len(split_annots[i]) == len(dst_indexes), \
                   "split/dst_indexes mismatch: %d %s %s" % (
                i, split_annots[i], len(dst_indexes))
        for i in range(len(dst_indexes)):
            h = MATRIX.headers_h[dst_indexes[i]]
            dst_annots = MATRIX.header2annots[h]
            assert len(split_annots) == len(dst_annots)
            for j in range(len(split_annots)):
                # change in place
                dst_annots[j] = split_annots[j][i]
    return MATRIX


def vcf_calc_vaf(MATRIX, calc_vaf):
    # List of: <ref index>,<alt index>,<vaf index>
    if not calc_vaf:
        return MATRIX

    jobs = []
    for x in calc_vaf:
        x = x.split(",")
        assert len(x) == 3, \
               "format should be: <ref index>,<alt index>,<vaf index>"
        x1, x2, x3 = x
        x1 = parse_indexes(MATRIX, x1)
        x2 = parse_indexes(MATRIX, x2)
        x3 = parse_indexes(MATRIX, x3)
        assert len(x1) == 1, x1
        assert len(x2) == 1, x2
        assert len(x3) == 1, x3
        ref_index, alt_index, vaf_index = x1[0], x2[0], x3[0]
        x = ref_index, alt_index, vaf_index
        jobs.append(x)
        
    MATRIX = MATRIX.copy()
    for x in jobs:
        ref_index, alt_index, vaf_index = x

        ref_h = MATRIX.headers_h[ref_index]
        alt_h = MATRIX.headers_h[alt_index]
        vaf_h = MATRIX.headers_h[vaf_index]
        ref_annots = MATRIX.header2annots[ref_h]
        alt_annots = MATRIX.header2annots[alt_h]
        vaf_annots = MATRIX.header2annots[vaf_h]

        # Change MATRIX in place.
        for i in range(len(ref_annots)):
            r = ref_annots[i].strip()
            a = alt_annots[i].strip()
            if not r or not a:
                continue
            r, a = int(r), int(a)
            total = r+a
            if not total:
                continue
            vaf_annots[i] = a / float(total)
    return MATRIX


def subtract_two_bed_lists(MATRIX, subtract_two_bed_lists):
    # Format: <annot 1>,<annot 2>,<dest>.  Each are 1-based
    # indexes.  <annot 1> is comma-separated list of numbers.  May end
    # in an extra comma.
    if not subtract_two_bed_lists:
        return MATRIX

    # same as calcBlocksizes, but order of annots is reversed.  Should
    # we keep both?

    x = subtract_two_bed_lists.split(",")
    assert len(x) == 3, "format should be: <annot1>,<annot2>,<dest>"
    i_1, i_2, i_dest = x
    i_1, i_2, i_dest = int(i_1), int(i_2), int(i_dest)
    # Convert to 0-based index.
    i_1, i_2, i_dest = i_1-1, i_2-1, i_dest-1
    assert i_1 >= 0 and i_1 < len(MATRIX.headers)
    assert i_2 >= 0 and i_2 < len(MATRIX.headers)
    assert i_dest >= 0 and i_dest < len(MATRIX.headers)

    MATRIX = MATRIX.copy()
    h_1 = MATRIX.headers_h[i_1]
    h_2 = MATRIX.headers_h[i_2]
    h_dest = MATRIX.headers_h[i_dest]

    annots_1 = MATRIX.header2annots[h_1]
    annots_2 = MATRIX.header2annots[h_2]
    assert len(annots_1) == len(annots_2)
    annots_dest = [""] * len(annots_1)
    for i in range(len(annots_1)):
        a1 = annots_1[i]
        a2 = annots_2[i]
        if not a1.strip() or not a2.strip():
            continue
        ends_with_comma = False
        a1 = a1.split(",")
        a2 = a2.split(",")
        assert len(a1) == len(a2), "Unequal lengths"
        if a1[-1] == "" or a2[-1] == "":
            ends_with_comma = True
        if ends_with_comma:
            assert a1[-1] == ""
            assert a2[-1] == ""
            a1 = a1[:-1]
            a2 = a2[:-1]
        a1 = [int(x) for x in a1]
        a2 = [int(x) for x in a2]
        d = [(a1[j]-a2[j]) for j in range(len(a1))]
        d = ",".join(map(str, d))
        if ends_with_comma:
            d = d + ","
        annots_dest[i] = d
    MATRIX.header2annots[h_dest] = annots_dest
    return MATRIX

                            
def subtract_value_from_bed_list(MATRIX, subtract_value_from_bed_list):
    # Format: <annot 1>,<annot 2>,<dest>.  Each are 1-based
    # indexes.  <annot 1> is comma-separated list of numbers.  May end
    # in an extra comma.  <annot 2> is single value.
    if not subtract_value_from_bed_list:
        return MATRIX

    x = subtract_value_from_bed_list.split(",")
    assert len(x) == 3, "format should be: <annot1>,<annot2>,<dest>"
    i_1, i_2, i_dest = x
    i_1, i_2, i_dest = int(i_1), int(i_2), int(i_dest)
    # Convert to 0-based index.
    i_1, i_2, i_dest = i_1-1, i_2-1, i_dest-1
    assert i_1 >= 0 and i_1 < len(MATRIX.headers)
    assert i_2 >= 0 and i_2 < len(MATRIX.headers)
    assert i_dest >= 0 and i_dest < len(MATRIX.headers)

    MATRIX = MATRIX.copy()
    h_1 = MATRIX.headers_h[i_1]
    h_2 = MATRIX.headers_h[i_2]
    h_dest = MATRIX.headers_h[i_dest]

    annots_1 = MATRIX.header2annots[h_1]
    annots_2 = MATRIX.header2annots[h_2]
    assert len(annots_1) == len(annots_2)
    annots_dest = [""] * len(annots_1)
    for i in range(len(annots_1)):
        a1 = annots_1[i]
        a2 = annots_2[i]
        if not a1.strip() or not a2.strip():
            continue
        ends_with_comma = False
        a1 = a1.split(",")
        if a1[-1] == "":
            ends_with_comma = True
            a1 = a1[:-1]
        a1 = [int(x) for x in a1]
        a2 = int(a2)
        d = [x-a2 for x in a1]
        d = ",".join(map(str, d))
        if ends_with_comma:
            d = d + ","
        annots_dest[i] = d
    MATRIX.header2annots[h_dest] = annots_dest
    return MATRIX


## def calc_blockSizes(MATRIX, calc_blockSizes):
##     # Format: <annot 1>,<annot 2>,<dest>.  Each are 1-based
##     # indexes.  <annot 1> is comma-separated list of numbers.  May end
##     # in an extra comma.
##     if not calc_blockSizes:
##         return MATRIX

##     x = calc_blockSizes.split(",")
##     assert len(x) == 3, "format should be: <annot1>,<annot2>,<dest>"
##     i_1, i_2, i_dest = x
##     i_1, i_2, i_dest = int(i_1), int(i_2), int(i_dest)
##     # Convert to 0-based index.
##     i_1, i_2, i_dest = i_1-1, i_2-1, i_dest-1
##     assert i_1 >= 0 and i_1 < len(MATRIX.headers)
##     assert i_2 >= 0 and i_2 < len(MATRIX.headers)
##     assert i_dest >= 0 and i_dest < len(MATRIX.headers)

##     MATRIX = MATRIX.copy()
##     h_1 = MATRIX.headers_h[i_1]
##     h_2 = MATRIX.headers_h[i_2]
##     h_dest = MATRIX.headers_h[i_dest]

##     annots_1 = MATRIX.header2annots[h_1]
##     annots_2 = MATRIX.header2annots[h_2]
##     assert len(annots_1) == len(annots_2)
##     annots_dest = [""] * len(annots_1)
##     for i in range(len(annots_1)):
##         a1 = annots_1[i]
##         a2 = annots_2[i]
##         if not a1.strip() or not a2.strip():
##             continue
##         ends_with_comma = False
##         a1 = a1.split(",")
##         a2 = a2.split(",")
##         if a1[-1] == "" or a2[-1] == "":
##             ends_with_comma = True
##         if ends_with_comma:
##             assert a1[-1] == ""
##             assert a2[-1] == ""
##             a1 = a1[:-1]
##             a2 = a2[:-1]
##         a1 = [int(x) for x in a1]
##         a2 = [int(x) for x in a2]
##         d = [(a2[j]-a1[j]) for j in range(len(a1))]
##         d = ",".join(map(str, d))
##         if ends_with_comma:
##             d = d + ","
##         annots_dest[i] = d
##     MATRIX.header2annots[h_dest] = annots_dest
##     return MATRIX

def _int_or_float(x):
    EPS = 1E-10
    x1 = float(x)
    try:
        x2 = int(x)
    except ValueError, x:
        return x1
    x = x1
    if (x1-x2) < EPS:
        x = x2
    return x

FILENAME = None  # for debugging

def main():
    global FILENAME
    import sys
    import argparse
    from genomicode import AnnotationMatrix

    parser = argparse.ArgumentParser(
        description="Perform operations on an annotation file.")
    parser.add_argument("filename", nargs=1, help="Annotation file.")
    parser.add_argument(
        "--read_as_csv", action="store_true",
        help="Read as a CSV file.")
    parser.add_argument(
        "--ignore_lines_startswith",
        help="Ignore lines that starts with this string.  "
        'E.g. --ignore_lines_starswith "##" will ignore headers in VCF files.')

    group = parser.add_argument_group(title="Matrix operations")
    group.add_argument(
        "--indexes", "--cut", dest="indexes", default=[], action="append",
        help="Select only these indexes from the file e.g. 1-5,8 "
        "(1-based, inclusive).  (MULTI)")
    group.add_argument(
        "--select_cols_str", default=[], action="append",
        help="Select the columns whose header contains matches this string.  "
        "(MULTI)")
    group.add_argument(
        "--select_cols_substr", default=[], action="append",
        help="Select the columns whose header contains this substring.  "
        "(MULTI)")
    group.add_argument(
        "--add_column", default=[], action="append",
        help="Add one or more columns.  "
        "Format: <index>,<header>,<default value>.  The column will be "
        "added before <index> (1-based).  If <index> is 1, this will be "
        'the new first column.  If <index> is "END", this will be '
        "the last column.  (MULTI)")
    group.add_argument(
        "--copy_column", default=[], action="append",
        help="Copy a column.  Format: <index>,<new_header>.  (MULTI)")
    
    group = parser.add_argument_group(title="Changing headers")
    group.add_argument(
        "--add_header_line", default=[], action="append",
        help="Add a header line to a file with no headers.  "
        "Format: <header1>[,<header2>...].  (MULTI)")
    group.add_argument(
        "--fill_empty_headers", action="store_true",
        help="If the header line contains some blanks, fill them in with "
        "defaults.")
    group.add_argument(
        "--remove_header_line", action="store_true",
        help="Remove the header line from the file.")
    group.add_argument(
        "--reorder_headers_alphabetical", action="store_true",
        help="Change the order of the headers.")
    group.add_argument(
        "--upper_headers", action="store_true",
        help="Make headers upper case.")
    group.add_argument(
        "--lower_headers", action="store_true",
        help="Make headers lower case.")
    group.add_argument(
        "--hash_headers", action="store_true",
        help="Hash the names of the headers.")
    group.add_argument(
        "--remove_duplicate_headers", action="store_true",
        help="If a matrix contains columns with the same header, "
        "keep only the first column.")
    group.add_argument(
        "--rename_duplicate_headers", action="store_true",
        help="Make all the headers unique.")
    group.add_argument(
        "--rename_header", default=[], action="append",
        help="Rename a header.  Format: <from>,<to>.  "
        "<from> will be replaced with <to>.  (MULTI)")
    group.add_argument(
        "--rename_header_i", default=[], action="append",
        help="Rename a header.  Format: <index>,<to>.  "
        "<index> is a 1-based column index.  (MULTI)")
    group.add_argument(
        "--append_to_headers", default=[], action="append",
        help="Append text to one or more headers.  "
        "Format: <indexes>;<text_to_append>.  (MULTI)")
    group.add_argument(
        "--prepend_to_headers", default=[], action="append",
        help="Prepend text to one or more headers.  "
        "Format: <indexes>;<text_to_prepend>.  (MULTI)")
    group.add_argument(
        "--replace_header", default=[], action="append",
        help="Replace a (sub)string with another in all headers.  "
        "Format: <from>,<to>.  <from> will be replaced with <to>.  (MULTI)")
    group.add_argument(
        "--replace_header_re", default=[], action="append",
        help="Like replace_header, but <from> can be a regular expression.  "
        "Format: <from>,<to>.  <from> will be replaced with <to>.  (MULTI)")
    
    group = parser.add_argument_group(title="Changing Annotations")
    group.add_argument(
        "--strip_all_annots", action="store_true",
        help="Get rid of spaces around each of the annotations.")
    group.add_argument(
        "--upper_annots", 
        help="Convert annotations to upper case.  Format: 1-based indexes.")
    group.add_argument(
        "--lower_annots", 
        help="Convert annotations to lower case.  Format: 1-based indexes.")
    group.add_argument(
        "--set_value_if_empty", default=[], action="append",
        help="If the column is empty, set with this value.  "
        "Format: <index 1-based>,<value>.  (MULTI)")
    group.add_argument(
        "--copy_value_if_empty", default=[], action="append",
        help="If the dest column is empty, copy the value from the src "
        "columns.  "
        "Format: <dest col>,<src col 1>[, <src col 2>...].  Columns "
        "are given as 1-based indexes.  (MULTI)")
    group.add_argument(
        "--copy_value_if_empty_header", default=[], action="append",
        help="Fill empty annotations with values from other columns "
        "with this header.  Gets the value from the left-most non-empty "
        "column with the same header.  "
        "Format: <dest header>,<src header 1>[, <src header 2>...].  (MULTI)")
    group.add_argument(
        "--copy_value_if_empty_same_header", default=[], action="append",
        help="Fill empty annotations with values from other columns "
        "that share this header.  Gets the value from the left-most non-empty "
        "column with the same header.  (MULTI)")
    group.add_argument(
        "--copy_value_if_empty_same_header_all", action="store_true",
        help="Fill empty annotations with values from other columns "
        "that share the same header.  Do for all columns that share the same "
        "header.")
    group.add_argument(
        "--rename_annot", default=[], action="append",
        help="Replace one whole annotation (not a substring) with another.  "
        "Format: <indexes>;<src>;<dst>.  (MULTI)")
    group.add_argument(
        "--replace_annot", default=[], action="append",
        help="Replace a substring of an annotation with another substring.  "
        "Format: <indexes>;<src>;<dst>.  (MULTI)")
    group.add_argument(
        "--prepend_to_annots", default=[], action="append",
        help="Prepend text to the values in one or more columns.  "
        "Format: <indexes>;<text_to_prepend>.  (MULTI)")
    group.add_argument(
        "--apply_re_to_annots", default=[], action="append",
        help="Apply a regular expression to annots.  "
        "Format: <indexes>;<regular expression>.  (MULTI)")
    group.add_argument(
        "--merge_annots", default=[], action="append",
        help="Merge a multiple annotations into one string.  "
        "Format: <src indexes>;<dst index>;<merge char>.  (MULTI)")
    group.add_argument(
        "--merge_annots_to_new_col", default=[], action="append",
        help="Merge a multiple annotations into one string.  "
        "Format: <src indexes>;<dst name>;<merge char>.  (MULTI)")
    group.add_argument(
        "--split_annots", default=[], action="append",
        help="Split an annotation across columns.  "
        "Format: <src index>;<dst indexes>;<split char>.  "
        "There should be at least one dst index for each item split.  (MULTI)")
    group.add_argument(
        "--split_chr_start_end", default=[], action="append",
        help='Split a chromosome location string (e.g. "chr1:320117-320142") '
        "into separate colummns: <chrom> <start> <end>.  "
        "Format: <header>.  <header> may be the name of the header, "
        "or a 1-based index.  (MULTI)")
    
    group = parser.add_argument_group(title="Mathematical Operations")
    group.add_argument(
        "--flip01",
        help="Flip 0's to 1's and 1's to 0's.  "
        "Format: indexes of columns to flip.")
    group.add_argument(
        "--all_same",
        help="Sets a 0 or 1 depending on whether the values in <indexes> "
        "are all the same.  "
        "Format: <indexes>;<index dest>.  All indexes should be 1-based.")
    group.add_argument(
        "--min_annots",
        help="Calculate the minimum value across a set of annotations.  "
        "Format: <indexes>;<index dest>.  All indexes should be 1-based.")
    group.add_argument(
        "--max_annots",
        help="Calculate the maximum value across a set of annotations.  "
        "Format: <indexes>;<index dest>.  All indexes should be 1-based.")
    group.add_argument(
        "--add_to", default=[], action="append",
        help="Add a number to a column.  "
        "Format: <header>,<number>.  "
        "Header can be the name of the header or a 1-based index.  (MULTI)")
    group.add_argument(
        "--multiply_by", default=[], action="append",
        help="Multiply a column by a number.  "
        "Format: <index>,<number>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--log_base", default=[], action="append",
        help="Log a column with a specific base.  "
        "Format: <index>,<base>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--neg_log_base", default=[], action="append",
        help="Log a column with a specific base and multiply by -1.  "
        "Format: <index>,<base>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--add_two_annots", default=[], action="append",
        help="Add column 1 to column 2 and save to a third column.  "
        "Format: <index 1>,<index 2>,<index dest>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--subtract_two_annots", default=[], action="append",
        help="Subtract column 2 from column 1 and save to a third column.  "
        "Format: <index 1>,<index 2>,<index dest>.  "
        "<index dest> = <index 1> - <index 2>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--divide_two_annots", default=[], action="append",
        help="Divide one column by another and save to a third column.  "
        "Format: <index numerator>,<index denominator>,<index dest>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--divide_many_annots", default=[], action="append",
        help="Divide a list of columns (in place) by another.  "
        "Format: <indexes numerator>;<index denominator>.  "
        "All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--average_same_header", action="store_true",
        help="Average the annotations that have the same header.")
    group.add_argument(
        "--round", default=[], action="append",
        help="Round the values of a column to integers.  "
        "Format: <index>.  All indexes should be 1-based.  (MULTI)")
    group.add_argument(
        "--convert_percent_to_decimal", default=[], action="append",
        help='Remove "%%" (if necessary) and divide by 100.  '
        "Format: <index>.  All indexes should be 1-based.  (MULTI)")
    

    group = parser.add_argument_group(title="VCF files")
    group.add_argument(
        "--vcf_standardize", 
        help="Take a VCF file (from IACS, Platypus, or GATK) and put into "
        "a standard format.  "
        "Format:<info_header>,<format_header>[,<genotype_header>].  "
        "<genotype_header> is a list of optional headers for the genotype "
        "information.  If not given, will use the last-most columns in "
        "the file.")
    group.add_argument(
        "--vcf_remove_bad_coords", action="store_true",
        help="Somve VCF files contain bad start or end positions, "
        "e.g. 8e+07.  Maybe been through Excel?  Remove them.")
    group.add_argument(
        "--vcf_remove_multicalls", action="store_true",
        help="Take a VCF file in standard format and make sure there is "
        "only one call per variant.")
    group.add_argument(
        "--vcf_extract_format_values", 
        help="Take a VCF file and extract the values from the corresponding "
        "format column.  "
        "Format: <format header>,<values header>,<value>[,value].  "
        "Example: FORMAT,Sample1,RD,AD,DP,FREQ")
    group.add_argument(
        "--vcf_extract_info_values", 
        help="Take a VCF file and extract the values from the INFO "
        "column.  Creates new columns."
        "Format: <INFO header>,<value>[,value].  "
        "Example: INFO,TC,TR")
    group.add_argument(
        "--vcf_split_AD", default=[], action="append",
        help="Split the AD value across columns.  "
        "Format: <src index>;<dst indexes>.  "
        "There should be at least one dst index for each item split.")
    group.add_argument(
        "--vcf_calc_vaf", default=[], action="append",
        help="Calculate the variant allele frequency.  "
        "Format: <ref index>,<alt index>,<vaf index>.  (MULTI)")


    group = parser.add_argument_group(title="Application-Specific Stuff")
    ## group.add_argument(
    ##     "--calc_blockSizes", 
    ##     help="For BED files, calculate blockSizes from blockStarts and "
    ##     "blockEnds.  "
    ##     "Format: <blockStarts index>,<blockEnds index>,<index dest>.  "
    ##     "All indexes should be 1-based.")
    group.add_argument(
        "--subtract_two_bed_lists", 
        help="For BED files, subtract column 2 (comma-separated values) "
        "from column 1 (comma-separated values) and save to a third "
        "column (comma-separated values).  "
        "Format: <index 1>,<index 2>,<index dest>.  "
        "<index dest> = <index 1> - <index 2>.  "
        "All indexes should be 1-based.")
    group.add_argument(
        "--subtract_value_from_bed_list", 
        help="For BED files, subtract column 2 (one value) "
        "from column 1 (comma-separated values) and save to a third "
        "column (comma-separated values).  "
        "Format: <index 1>,<index 2>,<index dest>.  "
        "<index dest> = <index 1> - <index 2>.  "
        "All indexes should be 1-based.")


    args = parser.parse_args()
    assert len(args.filename) == 1
    FILENAME = args.filename[0]

    # Do operations that do not take a matrix.
    if args.add_header_line:
        MATRIX = add_header_line(
            args.filename[0], args.add_header_line, args.read_as_csv)
    elif args.remove_header_line:
        remove_header_line(args.filename[0], args.read_as_csv)
        sys.exit(0)
    else:
        # Read the matrix.
        MATRIX = AnnotationMatrix.read(
            args.filename[0], args.read_as_csv, args.ignore_lines_startswith)

    # Perform operations.
    MATRIX = indexes_matrix(MATRIX, args.indexes)
    MATRIX = select_cols_str(MATRIX, args.select_cols_str)
    MATRIX = select_cols_substr(MATRIX, args.select_cols_substr)
    MATRIX = add_column(MATRIX, args.add_column)
    MATRIX = copy_column(MATRIX, args.copy_column)

    # Changing the headers.
    MATRIX = fill_empty_headers(MATRIX, args.fill_empty_headers)
    MATRIX = reorder_headers_alphabetical(
        MATRIX, args.reorder_headers_alphabetical)
    MATRIX = upper_headers(MATRIX, args.upper_headers)
    MATRIX = lower_headers(MATRIX, args.lower_headers)
    MATRIX = hash_headers(MATRIX, args.hash_headers)
    MATRIX = remove_duplicate_headers(MATRIX, args.remove_duplicate_headers)
    MATRIX = rename_duplicate_headers(MATRIX, args.rename_duplicate_headers)
    MATRIX = rename_header(MATRIX, args.rename_header)
    MATRIX = rename_header_i(MATRIX, args.rename_header_i)
    MATRIX = replace_header(MATRIX, args.replace_header)
    MATRIX = replace_header_re(MATRIX, args.replace_header_re)
    MATRIX = append_to_headers(MATRIX, args.append_to_headers)
    MATRIX = prepend_to_headers(MATRIX, args.prepend_to_headers)

    # Changing the values.
    MATRIX = strip_all_annots(MATRIX, args.strip_all_annots)
    MATRIX = upper_annots(MATRIX, args.upper_annots)
    MATRIX = lower_annots(MATRIX, args.lower_annots)
    MATRIX = set_value_if_empty(MATRIX, args.set_value_if_empty)
    MATRIX = copy_value_if_empty(MATRIX, args.copy_value_if_empty)
    MATRIX = copy_value_if_empty_header(
        MATRIX, args.copy_value_if_empty_header)
    MATRIX = copy_value_if_empty_same_header(
        MATRIX, args.copy_value_if_empty_same_header)
    MATRIX = copy_value_if_empty_same_header_all(
        MATRIX, args.copy_value_if_empty_same_header_all)
    MATRIX = replace_annot(MATRIX, args.replace_annot)
    MATRIX = replace_whole_annot(MATRIX, args.rename_annot)
    MATRIX = prepend_to_annots(MATRIX, args.prepend_to_annots)
    MATRIX = apply_re_to_annots(MATRIX, args.apply_re_to_annots)
    MATRIX = merge_annots(MATRIX, args.merge_annots)
    MATRIX = merge_annots_to_new_col(MATRIX, args.merge_annots_to_new_col)
    MATRIX = split_annots(MATRIX, args.split_annots)
    MATRIX = split_chr_start_end(MATRIX, args.split_chr_start_end)

    # Math operations.
    MATRIX = flip01_matrix(MATRIX, args.flip01)
    MATRIX = all_same(MATRIX, args.all_same)
    MATRIX = min_annots(MATRIX, args.min_annots)
    MATRIX = max_annots(MATRIX, args.max_annots)
    MATRIX = log_base(MATRIX, args.log_base)
    MATRIX = neg_log_base(MATRIX, args.neg_log_base)
    MATRIX = add_to(MATRIX, args.add_to)
    MATRIX = multiply_by(MATRIX, args.multiply_by)
    MATRIX = add_two_annots(MATRIX, args.add_two_annots)
    MATRIX = subtract_two_annots(MATRIX, args.subtract_two_annots)
    MATRIX = divide_two_annots(MATRIX, args.divide_two_annots)
    MATRIX = divide_many_annots(MATRIX, args.divide_many_annots)
    MATRIX = average_same_header(MATRIX, args.average_same_header)
    MATRIX = round_annots(MATRIX, args.round)
    MATRIX = convert_percent_to_decimal(
        MATRIX, args.convert_percent_to_decimal)

    # VCF
    MATRIX = vcf_standardize(MATRIX, args.vcf_standardize)
    MATRIX = vcf_remove_bad_coords(MATRIX, args.vcf_remove_bad_coords)
    MATRIX = vcf_remove_multicalls(MATRIX, args.vcf_remove_multicalls)
    MATRIX = vcf_extract_format_values(MATRIX, args.vcf_extract_format_values)
    MATRIX = vcf_extract_info_values(MATRIX, args.vcf_extract_info_values)
    MATRIX = vcf_split_AD(MATRIX, args.vcf_split_AD)
    MATRIX = vcf_calc_vaf(MATRIX, args.vcf_calc_vaf)

    # Application-specific stuff
    #MATRIX = calc_blockSizes(MATRIX, args.calc_blockSizes)
    MATRIX = subtract_two_bed_lists(MATRIX, args.subtract_two_bed_lists)
    MATRIX = subtract_value_from_bed_list(
        MATRIX, args.subtract_value_from_bed_list)

    # Write the matrix back out.
    AnnotationMatrix.write(sys.stdout, MATRIX)

if __name__ == '__main__':
    main()
