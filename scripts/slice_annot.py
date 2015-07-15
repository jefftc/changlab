#!/usr/bin/env python

# Functions:
# read_annot
# write_annot
#
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
# strip_all_annots
# upper_annots
# lower_annots
# replace_annots
# replace_whole_annot
# prepend_to_annots
# apply_re_to_annots
# divide_two_annots

import os
import sys


# TODO: merge with aligh_matrices
class AnnotationMatrix:
    def __init__(self, headers, headers_h, header2annots):
        # headers is a list of the original headers.
        # headers_h are the headers, hashed to ensure uniqueness.
        # headers2annots is a dictionary of hashed headers to the list
        # of annotations.
        assert headers
        assert headers_h
        assert len(headers) == len(headers_h)
        assert sorted(headers_h) == sorted(header2annots)
        for x in headers_h[1:]:
            assert len(header2annots[x]) == len(header2annots[headers_h[0]])
        self.headers = headers[:]
        self.headers_h = headers_h[:]
        self.header2annots = header2annots.copy()
    def copy(self):
        return AnnotationMatrix(
            self.headers, self.headers_h, self.header2annots)


def _hash_headers_unique(headers):
    # Make sure the headers in all_headers is unique.
    header2I = {}  # header -> list of indexes
    for i, header in enumerate(headers):
        if header not in header2I:
            header2I[header] = []
        header2I[header].append(i)

    nodup = headers[:]
    for (header, I) in header2I.iteritems():
        if len(I) < 2:
            continue
        for i in range(len(I)):
            nodup[I[i]] = "%s_%d" % (header, i+1)
    return nodup


def read_annot(filename, is_csv):
    # Everything are strings.  No numeric conversion.
    import re
    from genomicode import genesetlib

    delimiter = "\t"
    if is_csv:
        delimiter = ","

    all_headers, all_annots = [], []
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True,
        delimiter=delimiter):
        name, description, annots = x

        # Hack: Some files contain special characters, which mess up
        # alignment. Fix this here.
        # na\xc3\xafve-WIBR3.5 hESC
        # na\xe2\x80\x9a\xc3\xa0\xc3\xb6\xe2\x88\x9a\xc3\xb2ve-C1.2 hiPSC
        annots = [re.sub("na\\W+ve", "naive", x) for x in annots]

        all_headers.append(name)
        all_annots.append(annots)

    headers_h = _hash_headers_unique(all_headers)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    return AnnotationMatrix(all_headers, headers_h, header2annots)


def write_annot(handle_or_file, annot_matrix):
    from genomicode import jmath
    matrix = []
    for i, header_h in enumerate(annot_matrix.headers_h):
        header = annot_matrix.headers[i]
        annots = annot_matrix.header2annots[header_h]
        x = [header] + annots
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)

    handle = handle_or_file
    if type(handle) is type(""):
        handle = open(handle, 'w')
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


def parse_indexes(MATRIX, indexes_str):
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

    # Remove duplicated indexes.  Need to preserve order.
    nodup = []
    for i in I:
        if i not in nodup:
            nodup.append(i)
    I = nodup

    return I


def indexes_matrix(MATRIX, indexes_list):
    # indexes is a list of strings indicating indexes.
    if not indexes_list:
        return MATRIX
    I = []
    for indexes in indexes_list:
        x = parse_indexes(MATRIX, indexes)
        I.extend(x)

    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
    headers = [MATRIX.headers[i] for i in I]
    headers_h = [MATRIX.headers_h[i] for i in I]
    header2annots = {}
    for header_h in headers_h:
        header2annots[header_h] = MATRIX.header2annots[header_h]
    return AnnotationMatrix(headers, headers_h, header2annots)


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


def rename_duplicate_headers(MATRIX, rename_dups):
    if not rename_dups:
        return MATRIX

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

    headers = nodup
    headers_h = _hash_headers_unique(headers)
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)
    

def rename_header(MATRIX, rename_list):
    # rename_list is list of strings in format of: <from>,<to>.
    if not rename_list:
        return MATRIX

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
    headers_h = _hash_headers_unique(headers)
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)
        

def rename_header_i(MATRIX, rename_list):
    # rename_list is list of strings in format of: <index>,<to>.
    if not rename_list:
        return MATRIX

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
    headers_h = _hash_headers_unique(headers)
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)


def replace_header(MATRIX, replace_list):
    # replace_list is list of strings in format of: <from>,<to>.
    if not replace_list:
        return MATRIX

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
    headers_h = _hash_headers_unique(headers)
    
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)
        

def replace_header_re(MATRIX, replace_list):
    # replace_list is list of strings in format of: <from re>,<to>.
    import re

    if not replace_list:
        return MATRIX

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
                x = x.replace(m.group(0), to_str)
            headers[i] = x
    headers_h = _hash_headers_unique(headers)
    
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)


def prepend_to_headers(MATRIX, prepend_to_headers):
    # prepend_to_headers is list of strings in format of: <indexes>;<prefix>.
    if not prepend_to_headers:
        return MATRIX

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
    headers_h = _hash_headers_unique(headers)
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    return AnnotationMatrix(headers, headers_h, header2annots)
    

def add_column(MATRIX, add_column):
    # add_column is list of strings in format of: <index>,<header>,<default>.
    if not add_column:
        return MATRIX

    num_annots = None
    for annots in MATRIX.header2annots.itervalues():
        if num_annots is None:
            num_annots = len(annots)
        assert num_annots == len(annots)

    add_all = []  # list of (0-based index, header, default_value)
    for x in add_column:
        x = x.split(",", 2)
        assert len(x) == 3
        index, header, default_value = x
        index = int(index) - 1
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
    headers_h = _hash_headers_unique(headers)

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

    return AnnotationMatrix(headers, headers_h, header2annots)


def copy_column(MATRIX, copy_column):
    # copy_column is in format of: <index>,<new_header>.
    if not copy_column:
        return MATRIX

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
    headers_h = _hash_headers_unique(headers)

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

    return AnnotationMatrix(headers, headers_h, header2annots)


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
    header2annots = {}
    for header_h, annots in MATRIX.header2annots.iteritems():
        annots = [x.strip() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix(MATRIX.headers, MATRIX.headers_h, header2annots)


def upper_annots(MATRIX, upper):
    if not upper:
        return MATRIX
    I = parse_indexes(MATRIX, upper)
    
    header2annots = MATRIX.header2annots.copy()
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
        header_h = MATRIX.headers_h[i]
        annots = MATRIX.header2annots[header_h]
        annots = [x.upper() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix(MATRIX.headers, MATRIX.headers_h, header2annots)


def lower_annots(MATRIX, lower):
    if not lower:
        return MATRIX
    I = parse_indexes(MATRIX, lower)
    
    header2annots = MATRIX.header2annots.copy()
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
        header_h = MATRIX.headers_h[i]
        annots = MATRIX.header2annots[header_h]
        annots = [x.lower() for x in annots]
        header2annots[header_h] = annots
    return AnnotationMatrix(MATRIX.headers, MATRIX.headers_h, header2annots)


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


def subtract_two_annots(MATRIX, subtract_annots):
    # Format: <annot 1>,<annot 2>,<dest>.  Each are 1-based
    # indexes.  dest is annot1 - annot2.
    
    if not subtract_annots:
        return MATRIX

    x = subtract_annots.split(",")
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
        a1, a2 = float(a1), float(a2)
        annots_dest[i] = a1 - a2
    MATRIX.header2annots[h_dest] = annots_dest
    
    return MATRIX


def divide_two_annots(MATRIX, divide_annots):
    # Format: <numerator>,<denominator>,<dest>.  Each are 1-based
    # indexes.
    
    if not divide_annots:
        return MATRIX

    x = divide_annots.split(",")
    assert len(x) == 3, "format should be: <num>,<den>,<dest>"
    i_num, i_den, i_dest = x
    i_num, i_den, i_dest = int(i_num), int(i_den), int(i_dest)
    # Convert to 0-based index.
    i_num, i_den, i_dest = i_num-1, i_den-1, i_dest-1
    assert i_num >= 0 and i_num < len(MATRIX.headers)
    assert i_den >= 0 and i_den < len(MATRIX.headers)
    assert i_dest >= 0 and i_dest < len(MATRIX.headers)

    MATRIX = MATRIX.copy()
    h_num = MATRIX.headers_h[i_num]
    h_den = MATRIX.headers_h[i_den]
    h_dest = MATRIX.headers_h[i_dest]
    annots_num = MATRIX.header2annots[h_num]
    annots_den = MATRIX.header2annots[h_den]
    assert len(annots_num) == len(annots_den)
    annots_dest = [""] * len(annots_num)
    for i in range(len(annots_num)):
        num = annots_num[i]
        den = annots_den[i]
        if not num.strip() or not den.strip():
            continue
        num = float(num)
        den = float(den)
        # No divide by 0.
        if abs(den) < 1E-50:
            continue
        annots_dest[i] = num / den
    MATRIX.header2annots[h_dest] = annots_dest
    
    return MATRIX


def main():
    import argparse
    import arrayio

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
        "--add_column", default=[], action="append",
        help="Add one or more columns.  "
        "Format: <index>,<header>,<default value>.  The column will be "
        "added before <index> (1-based).  If <index> is 1, this will be "
        "the new first column.  (MULTI)")
    group.add_argument(
        "--copy_column", 
        help="Copy a column.  Format: <index>,<new_header>.")
    
    group = parser.add_argument_group(title="Changing headers")
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
        "--flip01",
        help="Flip 0's to 1's and 1's to 0's.  "
        "Format: indexes of columns to flip.")
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
        "--subtract_two_annots", 
        help="Subtract column 2 from column 1 and save to a third column.  "
        "Format: <index 1>,<index 2>,<index dest>.  "
        "All indexes should be 1-based.")
    group.add_argument(
        "--divide_two_annots", 
        help="Divide one column by another and save to a third column.  "
        "Format: <index numerator>,<index denominator>,<index dest>.  "
        "All indexes should be 1-based.")

    args = parser.parse_args()
    assert len(args.filename) == 1

    # Read the matrix.
    MATRIX = read_annot(args.filename[0], args.read_as_csv)

    # Perform operations.
    MATRIX = indexes_matrix(MATRIX, args.indexes)
    MATRIX = add_column(MATRIX, args.add_column)
    MATRIX = copy_column(MATRIX, args.copy_column)

    # Changing the headers.
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
    MATRIX = flip01_matrix(MATRIX, args.flip01)
    MATRIX = copy_value_if_empty(MATRIX, args.copy_value_if_empty)
    MATRIX = copy_value_if_empty_header(
        MATRIX, args.copy_value_if_empty_header)
    MATRIX = replace_annot(MATRIX, args.replace_annot)
    MATRIX = replace_whole_annot(MATRIX, args.rename_annot)
    MATRIX = prepend_to_annots(MATRIX, args.prepend_to_annots)
    MATRIX = apply_re_to_annots(MATRIX, args.apply_re_to_annots)
    MATRIX = divide_two_annots(MATRIX, args.divide_two_annots)
    MATRIX = subtract_two_annots(MATRIX, args.subtract_two_annots)

    # Write the matrix back out.
    write_annot(sys.stdout, MATRIX)

if __name__ == '__main__':
    main()
