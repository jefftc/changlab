#!/usr/bin/env python

# Functions:
# parse_indexes
#
# indexes_matrix
# flip01_matrix
#
# rename_duplicate_headers
# rename_header
# rename_header_i
# replace_header
# replace_header_re
# prepend_to_headers
#
# add_column
# copy_column
# 
# copy_value_if_empty
# copy_value_if_empty_header
# copy_value_if_empty_same_header
# strip_all_annots
# upper_annots
# lower_annots
# replace_annots
# replace_whole_annot
# prepend_to_annots
# apply_re_to_annots
#
# add_two_annots
# subtract_two_annots
# divide_two_annots
# divide_many_annots
# average_same_header


def parse_indexes(MATRIX, indexes_str, allow_duplicates=False):
    # Takes 1-based indexes and returns a list of 0-based indexes.
    # 
    # Example inputs:
    # 5
    # 1,5,10
    # 1-99,215-300
    from genomicode import parselib

    max_index = len(MATRIX.headers)

    I = []
    for s, e in parselib.parse_ranges(indexes_str):
        assert s >= 1
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
    import StringIO
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
    assert len(x) == len(X), "Matrix has %d columns, but %d headers given." % (
        len(X), len(x))
    headers = x
    headers_h = AnnotationMatrix.uniquify_headers(headers)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, X):
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
    # copy_column is in format of: <index>,<new_header>.
    if not copy_column:
        return MATRIX
    from genomicode import AnnotationMatrix

    x = copy_column.split(",", 1)
    assert len(x) == 2
    index, new_header = x
    index = int(index)
    assert index >= 1 and index <= len(MATRIX.headers)
    index -= 1  # convert to 0-based

    # Since the hashed header names might change, keep track of the
    # indexes for each header.
    h_indexes = [("OLD", i) for i in range(len(MATRIX.headers))]
    assert index >= 0 and index < len(h_indexes)
    h_indexes.insert(index+1, ("NEW", 0))

    headers = []
    for (which_one, i) in h_indexes:
        if which_one == "OLD":
            headers.append(MATRIX.headers[i])
        elif which_one == "NEW":
            headers.append(new_header)
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
            annots = MATRIX.header2annots[MATRIX.headers_h[index]]
            new_header_h = headers_h[i_new]
            header2annots[new_header_h] = annots
        else:
            raise AssertionError

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
    for indexes in copy_indexes:
        for x in itertools.product(indexes, indexes):
            i_dst, i_src = x
            if i_dst == i_src:
                continue
            header_dst = MATRIX.headers_h[i_dst]
            header_src = MATRIX.headers_h[i_src]
            # Change the annotations in place.
            annots_dst = MATRIX.header2annots[header_dst]
            annots_src = MATRIX.header2annots[header_src]
            for i in range(len(annots_dst)):
                # If dst is empty, then copy from src.
                if not annots_dst[i].strip() and annots_src[i].strip():
                    annots_dst[i] = annots_src[i]
    return MATRIX


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
        for i in range(MATRIX.num_annots()):
            all_annots = []
            for s_i in src_indexes:
                h = MATRIX.headers_h[s_i]
                annots = MATRIX.header2annots[h]
                annot = annots[i]
                all_annots.append(annot)
            merged = merge_char.join(all_annots)
            h = MATRIX.headers_h[dst_index]
            annots = MATRIX.header2annots[h]
            annots[i] = merged
    return MATRIX


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
            a1, a2 = float(a1), float(a2)
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
    # Format: list of <numerator indexes>;<denominator index.  Each
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


def main():
    import sys
    import argparse
    from genomicode import AnnotationMatrix

    parser = argparse.ArgumentParser(
        description="Perform operations on an annotation file.")
    parser.add_argument("filename", nargs=1, help="Annotation file.")
    parser.add_argument(
        "--read_as_csv", action="store_true",
        help="Read as a CSV file.")

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
        "--copy_column", 
        help="Copy a column.  Format: <index>,<new_header>.")
    
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
        "--hash_headers", action="store_true",
        help="Hash the names of the headers.")
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

    # Do operations that do not take a matrix.
    if args.add_header_line:
        MATRIX = add_header_line(
            args.filename[0], args.add_header_line, args.read_as_csv)
    elif args.remove_header_line:
        remove_header_line(args.filename[0], args.read_as_csv)
        sys.exit(0)
    else:
        # Read the matrix.
        MATRIX = AnnotationMatrix.read(args.filename[0], args.read_as_csv)

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
    MATRIX = hash_headers(MATRIX, args.hash_headers)
    MATRIX = rename_duplicate_headers(MATRIX, args.rename_duplicate_headers)
    MATRIX = rename_header(MATRIX, args.rename_header)
    MATRIX = rename_header_i(MATRIX, args.rename_header_i)
    MATRIX = replace_header(MATRIX, args.replace_header)
    MATRIX = replace_header_re(MATRIX, args.replace_header_re)
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
    MATRIX = replace_annot(MATRIX, args.replace_annot)
    MATRIX = replace_whole_annot(MATRIX, args.rename_annot)
    MATRIX = prepend_to_annots(MATRIX, args.prepend_to_annots)
    MATRIX = apply_re_to_annots(MATRIX, args.apply_re_to_annots)
    MATRIX = merge_annots(MATRIX, args.merge_annots)

    # Math operations.
    MATRIX = flip01_matrix(MATRIX, args.flip01)
    MATRIX = all_same(MATRIX, args.all_same)
    MATRIX = min_annots(MATRIX, args.min_annots)
    MATRIX = max_annots(MATRIX, args.max_annots)
    MATRIX = log_base(MATRIX, args.log_base)
    MATRIX = multiply_by(MATRIX, args.multiply_by)
    MATRIX = add_two_annots(MATRIX, args.add_two_annots)
    MATRIX = subtract_two_annots(MATRIX, args.subtract_two_annots)
    MATRIX = divide_two_annots(MATRIX, args.divide_two_annots)
    MATRIX = divide_many_annots(MATRIX, args.divide_many_annots)
    MATRIX = average_same_header(MATRIX, args.average_same_header)

    # Application-specific stuff
    #MATRIX = calc_blockSizes(MATRIX, args.calc_blockSizes)
    MATRIX = subtract_two_bed_lists(MATRIX, args.subtract_two_bed_lists)
    MATRIX = subtract_value_from_bed_list(
        MATRIX, args.subtract_value_from_bed_list)

    # Write the matrix back out.
    AnnotationMatrix.write(sys.stdout, MATRIX)

if __name__ == '__main__':
    main()
