#!/usr/bin/env python

# Functions:
# read_matrices
# 
# parse_indexes
# parse_names
# parse_geneset
# _parse_file_gs
# _parse_file_annot
# _read_annot_file
#
# find_col_indexes
# find_col_genesets
# find_col_annotation
# relabel_col_ids
# remove_duplicate_cols
# remove_col_ids
# 
# find_row_indexes
# find_row_ids
# find_row_genesets
# find_row_annotation
# find_row_numeric_annotation
# find_row_mean_var
# 
# align_rows
# align_cols
# 
# add_row_annot
# remove_row_annot
# rename_row_annot
# move_row_annot
#
# _match_rownames_to_geneset       DEPRECATED
# _match_colnames_to_geneset       DEPRECATED
# _best_match_colnames_to_geneset  DEPRECATED
# _num_matches_to_geneset          DEPRECATED
# _align_geneset_to_matrix         DEPRECATED
# _align_matrix_to_geneset         DEPRECATED
# _find_col_header                 DEPRECATED
# 
# _intersect_indexes
# _dedup_indexes

import os, sys

def _clean(s):
    s = s.replace("\t", " ")
    s = s.strip()
    if s.startswith('"') and s.endswith('"'):
        s = s[1:-1]
    return s

def read_matrices(filenames, read_as_csv, remove_comments, clean_only):
    import csv
    import tempfile
    import arrayio
    from genomicode import filelib
    
    temp_files = []
    try:
        if read_as_csv or (remove_comments is not None):
            delimiter = "\t"
            if read_as_csv:
                delimiter = ","
            assert delimiter not in remove_comments
                
            # Make a temporary files for each matrix file.
            for filename in filenames:
                x, file = tempfile.mkstemp(dir="."); os.close(x)
                temp_files.append(file)
            assert len(filenames) == len(temp_files)

            # Clean up the files.
            for (infile, outfile) in zip(filenames, temp_files):
                if read_as_csv:
                    inhandle = csv.reader(filelib.openfh(infile))
                else:
                    inhandle = filelib.read_cols(infile, delimiter=delimiter)
                    
                outhandle = open(outfile, 'w')
                for cols in inhandle:
                    if not cols:
                        continue
                    if remove_comments and cols[0].startswith(remove_comments):
                        continue
                    # Clean up the data.
                    cols = [_clean(x) for x in cols]
                    
                    outhandle.write("\t".join(cols)+"\n")
                outhandle.close()

            # Use the cleaned files.
            filenames = temp_files

        if clean_only:
            if len(filenames) != 1:
                raise NotImplementedError, "Can not clean and merge."
            sys.stdout.write(open(filenames[0]).read())
            sys.exit(0)

        # Read the files.
        matrices = []
        for filename in filenames:
            fmt_module = arrayio.guess_format(filename)
            assert fmt_module,"I could not figure out the format of file: %s"%\
                filename
            x = fmt_module.read(filename)
            matrices.append(x)
    finally:
        for file in temp_files:
            if os.path.exists(file):
                os.unlink(file)
    
    return fmt_module, matrices

def parse_indexes(MATRIX, is_row, s):
    # Examples:
    # 5
    # 1,5,10
    # 1-99,215-300
    from genomicode import parselib

    max_index = MATRIX.nrow()
    if not is_row:
        max_index = MATRIX.ncol()
    
    I = []
    for s, e in parselib.parse_ranges(s):
        assert s >= 1
        s, e = s-1, min(e, max_index)
        I.extend(range(s, e))
    return I

def parse_names(MATRIX, is_row, s):
    # Examples:
    # E2F1,E2F3,PCNA
    names = s.split(",")
    params = { "row" : names }
    if not is_row:
        params = { "col" : names }
        
    I_row, I_col = MATRIX._index(**params)
    I = I_row
    if not is_row:
        I = I_col
    return I

def parse_geneset(MATRIX, is_row, geneset):
    # Return a list of indexes that match the desired gene sets.
    # If no genesets are specified, return the indexes that match any
    # of the genesets.
    from genomicode import genesetlib

    if not geneset:
        return []
    filename, genesets = _parse_file_gs(geneset)

    keywds = { "allow_tdf" : True }
    genes = genesetlib.read_genes(filename, *genesets, **keywds)
    params = { "row" : genes }
    if not is_row:
        params = { "col" : genes }
    I_row, I_col = MATRIX._index(**params)
    I = I_row
    if not is_row:
        I = I_col
    return I

def _parse_file_gs(geneset):
    # Parse a geneset specified by the user.  geneset is in the format
    # of <filename>[,<geneset>,<geneset>,...].  Return a tuple of
    # <filename>, list of <geneset> (or empty list).
    # XXX what happens if this is an empty list?
    x = geneset.split(",")
    assert len(x) >= 1
    filename, genesets = x[0], x[1:]
    return filename, genesets

def _parse_file_annot(annotation):
    # annotations is in the format:
    # <txt_file>,<header>,<value>[,<value,...]
    # Return a tuple of <filename>, <header>, list of <value>.
    x = annotation.split(",")
    assert len(x) >= 3
    filename, header, values = x[0], x[1], x[2:]
    return filename, header, values

def _read_annot_file(filename):
    # Return (header2annots, all_headers, all_annots).
    from genomicode import genesetlib
    
    assert os.path.exists(filename)
    header2annots = {}
    all_headers = []
    all_annots = []
    for x in genesetlib.read_tdf(filename):
        name, x, annots = x
        header2annots[name] = annots
        all_headers.append(name)
        all_annots.append(annots)
    assert all_headers
    for annots in all_annots[1:]:
        assert len(annots) == len(all_annots[0])
    return header2annots, all_headers, all_annots

def _parse_file_num_annot(annotation):
    # annotations is in the format:
    # <txt_file>,<header>,<value>[,<value,...]
    # <value> can optionally start with modifiers "<", "<=", ">=" or ">".
    # Return a tuple of <filename>, <header>, list of (<modifier>, <value>).
    # If no modifier is specified, then <modifier> is "=".  <value> is
    # a floating point number.
    x = annotation.split(",")
    assert len(x) >= 3
    filename, header, values = x[0], x[1], x[2:]
    for i in range(len(values)):
        value = values[i]
        modifier = "="
        if value[:2] in ["<=", ">="]:
            modifier = value[:2]
            value = value[2:]
        elif value[0] in "<>=":
            modifier = value[0]
            value = value[1:]
        value = float(value)
        values[i] = modifier, value
    return filename, header, values
    
def find_col_indexes(MATRIX, indexes):
    if not indexes:
        return None
    return parse_indexes(MATRIX, False, indexes)

def find_col_genesets(MATRIX, genesets):
    if not genesets:
        return None
    return parse_geneset(MATRIX, False, genesets)

def find_col_annotation(MATRIX, col_annotation):
    # Format: <txt_file>,<header>,<value>[,<value,...]
    from genomicode import matrixlib
    
    if not col_annotation:
        return None

    x = _parse_file_annot(col_annotation)
    filename, header, values = x

    # Read the annotations.
    x = _read_annot_file(filename)
    header2annots, all_headers, all_annots = x

    # Align the annotations to the matrix file.
    x = matrixlib.align_cols_to_many_annots(
        MATRIX, all_annots, hash=True, get_indexes=True)
    I_matrix, I_annots, index = x
    assert len(I_matrix) == len(I_annots)

    annots = header2annots[header]
    I = []
    for (im, ia) in zip(I_matrix, I_annots):
        if annots[ia] in values:
            I.append(im)
    #I = [im for (im, ia) in zip(I_matrix, I_annots) if annots[ia] in values]

##     # Search for a col_annot in the MATRIX that matches a column in
##     # the annotation file.
##     for annot_header in header_order:
##         matrix_header = _find_col_header(MATRIX, header2annots[annot_header])
##         if matrix_header is not None:
##             break
##     else:
##         raise AssertionError, "I could not align the annotation file to " + \
##               "the matrix."
##     annot = header2annots[annot_header]
##     matrix_annot = MATRIX.col_names(matrix_header)
##     h_annot = [hashlib.hash_R(x) for x in annot]
##     h_matrix_annot = [hashlib.hash_R(x) for x in matrix_annot]
##     I = jmath.match(h_matrix_annot, h_annot)
##     assert None not in I
##     for h in header2annots:
##         x = header2annots[h]
##         header2annots[h] = [x[i] for i in I]
##     # Make sure they're aligned correctly.
##     annot = header2annots[annot_header]
##     matrix_annot = MATRIX.col_names(matrix_header)
##     h_annot = [hashlib.hash_R(x) for x in annot]
##     h_matrix_annot = [hashlib.hash_R(x) for x in matrix_annot]
##     assert h_annot == h_matrix_annot

##     # Identify the lines that match the annotations.
##     assert header in header2annots, "I could not find the header: %s" % header
##     I = []
##     for i, annot in enumerate(header2annots[header]):
##         if annot in values:
##             I.append(i)
    return I

def relabel_col_ids(MATRIX, geneset):
    import arrayio
    from genomicode import genesetlib
    from genomicode import matrixlib

    if not geneset:
        return MATRIX
    filename, genesets = _parse_file_gs(geneset)
    assert len(genesets) == 1
    print filename
    print genesets
    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    for x in genesetlib.read_genesets(filename, allow_tdf=True):
        geneset, description, genes = x
        geneset2genes[geneset] = genes
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Find an alignment between the sample names and the genesets.
    x = matrixlib.align_cols_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True)
    I_matrix, I_geneset, index = x
    if len(I_matrix) != MATRIX.ncol():
        # No matches.  Try to diagnose.
        #gs, nm, m, mm = _best_match_colnames_to_geneset(
        #    MATRIX, all_genesets, geneset2genes, hash=True)
        if len(I_matrix) == 0:
            print "Matrixes doesn't match any gene sets."
        else:
            gs = all_genesets[index]
            nm = len(I_matrix)
            missing = []
            col_names = MATRIX.col_names(arrayio.COL_ID)
            missing = [col_names[i] for i in range(len(col_names))
                       if i not in I_matrix]
            print "Matrix best matches %s [%d:%d]." % (gs, nm, MATRIX.ncol())
            if len(missing) > 5:
                print "Missing (showing %d of %d):" % (5, len(missing))
            else:
                print "Missing:"
            for x_ in missing[:5]:
                print x_
        raise AssertionError, "I could not match the matrix to a geneset."

    # Add the new column names to the MATRIX.
    MATRIX_new = MATRIX.matrix()
    name = arrayio.COL_ID
    if name not in MATRIX_new._col_names:
        name = MATRIX_new._synonyms[name]
    assert name in MATRIX_new._col_names, "I can not find the sample names."
    gs = genesets[0]
    x = geneset2genes[gs]
    names = [x[i] for i in I_geneset]
    MATRIX_new._col_names[name] = names
    
    return MATRIX_new

def remove_duplicate_cols(MATRIX, filter_duplicate_cols):
    import arrayio
    
    if not filter_duplicate_cols:
        return MATRIX
    headers = MATRIX.col_names(arrayio.COL_ID)
    I = []
    seen = {}
    for i in range(len(headers)):
        if headers[i] in seen:
            continue
        seen[headers[i]] = 1
        I.append(i)
    x = MATRIX.matrix(None, I)
    return x

def remove_col_ids(MATRIX, remove_col_ids):
    import arrayio

    if remove_col_ids is None:
        return MATRIX
    x = remove_col_ids.split(",")
    col_ids = [x.strip() for x in x]
    names = MATRIX.col_names(arrayio.COL_ID)
    
    I = []
    not_found = {}.fromkeys(col_ids)
    for i, name in enumerate(names):
        if name in col_ids:
            if name in not_found:
                del not_found[name]
        else:
            I.append(i)
    assert not not_found, "I could not find: %s." % \
        ", ".join(sorted(not_found))
    x = MATRIX.matrix(None, I)
    return x
    
def find_row_indexes(MATRIX, indexes):
    if not indexes:
        return None
    return parse_indexes(MATRIX, True, indexes)

def find_row_ids(MATRIX, ids):
    # IDs is a list of gene IDs or names.  Combine them into a single string.
    ids = [x.strip() for x in ids]
    ids = ",".join(ids)
    if not ids:
        return None
    return parse_names(MATRIX, True, ids)

def find_row_genesets(MATRIX, genesets):
    if not genesets:
        return None
    return parse_geneset(MATRIX, True, genesets)

def find_row_annotation(MATRIX, row_annotation):
    # Format: <txt_file>,<header>,<value>[,<value,...]
    from genomicode import matrixlib
    
    if not row_annotation:
        return None

    x = _parse_file_annot(row_annotation)
    filename, header, values = x

    # Read the annotations.
    x = _read_annot_file(filename)
    header2annots, all_headers, all_annots = x

    # Align the annotations to the matrix file.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_annots, hash=True, get_indexes=True)
    I_matrix, I_annots, index = x
    assert len(I_matrix) == len(I_annots)

    annots = header2annots[header]
    I = []
    for (im, ia) in zip(I_matrix, I_annots):
        if annots[ia] in values:
            I.append(im)
    return I

def find_row_numeric_annotation(MATRIX, row_annotation):
    # Format: <txt_file>,<header>,<value>[,<value,...]
    from genomicode import matrixlib
    
    if not row_annotation:
        return None

    x = _parse_file_num_annot(row_annotation)
    filename, header, values = x

    # Read the annotations.
    x = _read_annot_file(filename)
    header2annots, all_headers, all_annots = x

    # Align the annotations to the matrix file.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_annots, hash=True, get_indexes=True)
    I_matrix, I_annots, index = x
    assert len(I_matrix) == len(I_annots)

    annots = header2annots[header]
    I = []
    for (im, ia) in zip(I_matrix, I_annots):
        x = float(annots[ia])
        match = True
        for (modifier, value) in values:
            if modifier == "=":
                match = (x == value)
            elif modifier == "<":
                match = (x < value)
            elif modifier == ">":
                match = (x > value)
            elif modifier == "<=":
                match = (x <= value)
            elif modifier == ">=":
                match = (x >= value)
            else:
                raise AssertionError, "Unknown modifier: %s" % modifier
        if match:
            I.append(im)
    return I

def find_row_mean_var(MATRIX, filter_mean, filter_var):
    from genomicode import pcalib
    if filter_mean is None and filter_var is None:
        return None
    if filter_mean is not None:
        filter_mean = float(filter_mean)
    if filter_var is not None:
        filter_var = float(filter_var)
    
    assert filter_mean is None or (filter_mean >= 0 and filter_mean <= 1)
    assert filter_var is None or (filter_var >= 0 and filter_var <= 1)

    nrow = MATRIX.nrow()

    num_genes_mean = num_genes_var = None
    if filter_mean is not None:
        # Calculate the number of genes to keep.
        num_genes_mean = int((1.0 - filter_mean) * nrow)
    if filter_var is not None:
        # Calculate the number of genes to keep.
        num_genes_var = int((1.0 - filter_var) * nrow)
    I = pcalib.select_genes_mv(MATRIX._X, num_genes_mean, num_genes_var)
    return I

def align_rows(MATRIX, align_row_file, ignore_missing_rows):
    import arrayio
    
    if not align_row_file:
        return None
    assert os.path.exists(align_row_file), \
           "File not found: %s" % align_row_file

    ALIGN = arrayio.read(align_row_file)
    # Try all the headers and see if we can find a hit.
    # BUG: what if there's no header?
    best_I = []
    for header in ALIGN.row_names():
        ids = ALIGN.row_names(header)
        I_row, I_col = MATRIX._index(row=ids, row_header=arrayio.ROW_ID)
        # Bug: cannot just check length.  If there are duplicate rows,
        # the longest one may not be the one that matches the most
        # unique rows.
        if len(I_row) > len(best_I):
            best_I = I_row
        #I = I_row
        if len(best_I) == len(ids):
            break
    I = best_I
    if not ignore_missing_rows and len(ids) != len(I):
        # Diagnose problem here.
        x = ALIGN.row_names(arrayio.ROW_ID)
        ids_A = {}.fromkeys(x)
        x = MATRIX.row_names(arrayio.ROW_ID)
        ids_M = {}.fromkeys(x)
        missing = []  # In the align file, but not in my file.
        for id_ in ids_A:
            if id_ not in ids_M:
                missing.append(id)
        if len(missing) < 10:
            for id_ in sorted(missing):
                print id_
        message = ("%d IDs from the align file are missing from the "
                   "matrix file." % len(missing))
        raise AssertionError, message
    return I

def align_cols(MATRIX, align_col_file, ignore_missing_cols):
    import arrayio

    if not align_col_file:
        return None
    assert os.path.exists(align_col_file), \
        "File not found: %s" % align_col_file

    headers = MATRIX.col_names(arrayio.COL_ID)
    ALIGN = arrayio.read(align_col_file)
    # Try all the headers and see if we can find a hit.
    # Bug: what if there are duplicates in MATRIX or ALIGN?
    best_matches, best_I, best_header = None, [], ""
    for header in ALIGN.col_names():
        ids = ALIGN.col_names(header)
        I_row, I_col = MATRIX._index(col=ids, col_header=arrayio.COL_ID)
        # Count the number of unique matches.
        num_matches = {}
        for i in I_col:
            num_matches[headers[i]] = 1
        num_matches = len(num_matches)
        if best_matches is None or num_matches > best_matches:
            best_matches, best_I, best_header = num_matches, I_col, header
        if num_matches == len(ids):
            break
    I = best_I
    #if ignore_duplicate_cols:
    #    seen = {}
    #    i = 0
    #    while i < len(I):
    #        if headers[I[i]] in seen:
    #            del I[i]
    #        else:
    #            seen[headers[I[i]]] = 1
    #            i += 1
    #for i in best_I:
    #    print i, MATRIX.col_names(arrayio.COL_ID)[i]
    #sys.exit(0)
    if not ignore_missing_cols and len(ids) != len(I):
        # Diagnose problem here.
        x = ALIGN.col_names(arrayio.COL_ID)
        ids_A = {}.fromkeys(x)
        x = MATRIX.col_names(arrayio.COL_ID)
        ids_M = {}.fromkeys(x)
        missing = []
        for id_ in ids_A:
            if id_ not in ids_M:
                missing.append(id_)
        if len(missing) < 10:
            for id_ in sorted(missing):
                print id_
        message = "I could not find %d IDs." % len(missing)
        raise AssertionError, message
    return I

def add_row_annot(MATRIX, row_annots):
    # row_annot should be in the format <gmx/gmt_file>[,<geneset>].
    from genomicode import genesetlib
    from genomicode import matrixlib

    if not row_annots:
        return MATRIX
    filename, genesets = _parse_file_gs(row_annots)
    
    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    num_genes = None
    for x in genesetlib.read_genesets(filename, allow_tdf=True):
        geneset, description, genes = x
        geneset2genes[geneset] = genes
        if num_genes is None:
            num_genes = len(genes)
        assert len(genes) == num_genes, "%s %d %d" % (
            geneset, len(genes), num_genes)
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Find an alignment between one of the matrix row_names and the
    # genesets.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True)
    I_matrix, I_geneset, x = x
    assert len(I_matrix), "I could not match the matrix to a geneset."
    #x = _match_rownames_to_geneset(MATRIX, all_genesets, geneset2genes)
    #assert x, "I could not match the matrix to a geneset."
    #I_matrix, I_geneset = x
    MATRIX_new = MATRIX.matrix(I_matrix, None)

    # Add the new annotations to the MATRIX.
    for gs in genesets:
        assert gs in geneset2genes, "Missing geneset: %s" % gs
        assert gs not in MATRIX.row_names(), "duplicate name: %s" % gs
        
        x = geneset2genes[gs]
        genes = [x[i] for i in I_geneset]

        assert len(genes) == MATRIX_new.nrow()
        MATRIX_new._row_order.append(gs)
        MATRIX_new._row_names[gs] = genes
        
    return MATRIX_new
        
def remove_row_annot(MATRIX, name):
    assert name in MATRIX.row_names(), "I could not find name: %s" % name
    MATRIX_clean = MATRIX.matrix()
    Mc = MATRIX_clean
    assert name in Mc._row_names
    del Mc._row_names[name]
    if Mc._row_order and name in Mc._row_order:
        i = Mc._row_order.index(name)
        Mc._row_order.pop(i)
    # Bug: Does not remove from the synonyms list.
    return MATRIX_clean

def rename_row_annot(MATRIX, row_annot):
    if row_annot is None:
        return MATRIX
    assert "," in row_annot
    x = row_annot.split(",", 1)
    old_name, new_name = x
    assert old_name in MATRIX.row_names(), \
        "I could not find name: %s" % old_name
    
    MATRIX_clean = MATRIX.matrix()
    Mc = MATRIX_clean
    assert old_name in Mc._row_names
    Mc._row_names[new_name] = Mc._row_names[old_name]
    del Mc._row_names[old_name]
    if Mc._row_order and old_name in Mc._row_order:
        i = Mc._row_order.index(old_name)
        Mc._row_order[i] = new_name
    if hasattr(Mc, "_synonyms"):
        synonyms = {}
        for (key, value) in Mc._synonyms.iteritems():
            if key == old_name:
                key = new_name
            if value == old_name:
                value = new_name
            synonyms[key] = value
        Mc.__dict__["_synonyms"] = synonyms
    return MATRIX_clean

def move_row_annot(MATRIX, move_row_annot):
    if move_row_annot is None:
        return MATRIX
    x = move_row_annot.split(",", 1)
    assert len(x) == 2
    old_index, new_index = x
    old_index = int(old_index)
    new_index = int(new_index)

    MATRIX_clean = MATRIX.matrix()
    assert MATRIX_clean._row_order
    order = MATRIX_clean._row_order[:]
    assert old_index >= 1 and old_index <= len(order)
    assert new_index >= 1 and new_index <= len(order)
    x = order.pop(old_index-1)
    order.insert(new_index-1, x)
    MATRIX_clean._row_order = order
    return MATRIX_clean


## def _match_rownames_to_geneset(MATRIX, all_genesets, geneset2genes):
##     # Return tuple of (I_matrix, I_geneset) or None if no match can be
##     # found.  Will find the largest match possible.

##     # Align every geneset to every row name in the matrix.
##     geneset_aligns = []  # list of (I_geneset, rowname, geneset)
##     matrix_aligns = []   # list of (I_matrix, rowname, geneset)
##     for name in MATRIX.row_names():
##         annots = MATRIX.row_names(name)
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             I_geneset = _align_geneset_to_matrix(annots, genes)
##             I_matrix = _align_matrix_to_geneset(annots, genes)
##             if I_geneset:
##                 x = I_geneset, name, gs
##                 geneset_aligns.append(x)
##             if I_matrix:
##                 x = I_matrix, name, gs
##                 matrix_aligns.append(x)
                
##     # First, try to find a geneset that matches the exactly matrix.
##     # Favor geneset_aligns over matrix_aligns to avoid changing the
##     # matrix.
##     for x in geneset_aligns:
##         I_geneset, rowname, geneset = x
##         I_matrix = range(MATRIX.nrow())
##         assert len(I_matrix) == len(I_geneset)
##         return I_matrix, I_geneset

##     # Otherwise, choose the match that generates the largest matrix.
##     I_matrix = None
##     for x in matrix_aligns:
##         I, rowname, geneset = x
##         if I_matrix is None or len(I) > len(I_matrix):
##             I_matrix = I
##     if I_matrix:
##         I_geneset = range(len(I_matrix))
##         return I_matrix, I_geneset

##     return None

## def _match_colnames_to_geneset(
##     MATRIX, all_genesets, geneset2genes, hash=False):
##     # Return I_geneset or None if no match can be found.  The sample
##     # names from MATRIX must match the names in the gene set exactly.
##     # MATRIX may be a subset of the gene set.
##     import arrayio
##     from genomicode import hashlib

##     annots = MATRIX.col_names(arrayio.COL_ID)
##     if hash:
##         annots = [hashlib.hash_R(x) for x in annots]
##         g2g = {}
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             genes = [hashlib.hash_R(x) for x in genes]
##             g2g[gs] = genes
##         geneset2genes = g2g
        
##     # Align every geneset to the col names in the matrix.
##     geneset_aligns = []  # list of (I_geneset, geneset)
##     for gs in all_genesets:
##         genes = geneset2genes[gs]
##         I_geneset = _align_geneset_to_matrix(annots, genes)
##         if I_geneset:
##             x = I_geneset, gs
##             geneset_aligns.append(x)
                
##     # Find a geneset that exactly matches the sample names in the matrix.
##     for x in geneset_aligns:
##         I_geneset, geneset = x
##         # Allow the geneset to be a superset of MATRIX.
##         if len(I_geneset) > MATRIX.ncol():
##             continue
##         return I_geneset

##     # Return None if not found.
##     return None

## def _best_match_colnames_to_geneset(
##     MATRIX, all_genesets, geneset2genes, hash=False):
##     # Return the (name, num_matches, list of matches, list of
##     # non-matches) for the geneset that best matches the colnames.
##     import arrayio
##     from genomicode import hashlib

##     annots = MATRIX.col_names(arrayio.COL_ID)
##     annots_orig = annots[:]
##     if hash:
##         annots = [hashlib.hash_R(x) for x in annots]
##         g2g = {}
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             genes = [hashlib.hash_R(x) for x in genes]
##             g2g[gs] = genes
##         geneset2genes = g2g

##     # Align every geneset to the col names in the matrix.
##     best_gs = best_num_matches = best_matches = best_mismatches = None
##     for gs in all_genesets:
##         genes = geneset2genes[gs]
##         x = _num_matches_to_geneset(annots_orig, annots, genes)
##         num_matches, matches, mismatches = x
##         if best_num_matches is None or num_matches > best_num_matches:
##             best_gs = gs
##             best_num_matches = num_matches
##             best_matches = matches
##             best_mismatches = mismatches
                
##     return best_gs, best_num_matches, best_matches, best_mismatches

## def _num_matches_to_geneset(matrix_annots, hashed_annots, geneset_genes):
##     # Count the number of annotations that match the geneset.  Return
##     # number of matches, list of matches, list of mismatches.
##     assert len(matrix_annots) == len(hashed_annots)
    
##     count = 0
##     matches, mismatches = [], []
##     for i, annot in enumerate(hashed_annots):
##         if annot in geneset_genes:
##             count += 1
##             matches.append(matrix_annots[i])
##         else:
##             mismatches.append(matrix_annots[i])
##     return count, matches, mismatches

## def _align_geneset_to_matrix(matrix_annots, geneset_genes):
##     # Return a list of the indexes required to align the geneset to
##     # the rows of the matrix.  If it cannot be completely aligned,
##     # then return None.
##     from genomicode import jmath

##     if len(matrix_annots) > len(geneset_genes):
##         return None
##     I = jmath.match(matrix_annots, geneset_genes)
##     if None in I:
##         return None
##     return I

## def _align_matrix_to_geneset(matrix_annots, geneset_genes):
##     # Return a list of the indexes required to align the matrix to the
##     # geneset.  If it cannot be completely aligned, then return None.
##     from genomicode import jmath

##     if len(geneset_genes) > len(matrix_annots):
##         return None
##     I = jmath.match(geneset_genes, matrix_annots)
##     if None in I:
##         return None
##     return I

## def _find_col_header(MATRIX, col_names):
##     # Find the header in the MATRIX whose annotations can be found in
##     # col_names.  MATRIX can be a subset of col_names.  If multiple
##     # headers match, then just return the first one.  If none match,
##     # return None.
##     from genomicode import jmath
##     from genomicode import hashlib

##     h_col_names = [hashlib.hash_R(x) for x in col_names]

##     for name in MATRIX.col_names():
##         matrix_col_names = MATRIX.col_names(name)
##         h_matrix_col_names = [hashlib.hash_R(x) for x in matrix_col_names]
##         #I = jmath.match(h_col_names, h_matrix_col_names)
##         #if None in I:  # missing something
##         #    continue
##         I = jmath.match(h_matrix_col_names, h_col_names)
##         missing = []
##         for i in range(len(I)):
##             if I[i] is None:
##                 missing.append(matrix_col_names[i])
##         #print name, len(missing), missing
##         if len(missing):
##             continue
##         return name
##     return None

def _intersect_indexes(*indexes):
    # None means take all indexes.
    indexes = [x for x in indexes if x is not None]
    if not indexes:
        return None

    # Only want the indexes that occur in them all.  Preserve order.
    I = indexes[0]
    for ind in indexes[1:]:
        I = [i for i in I if i in ind]
    return I

def _dedup_indexes(I):
    # Get rid of duplicate indexes, preserving the original order of
    # the indexes.
    I = I[:]
    
    i = 0
    seen = {}
    while i < len(I):
        if I[i] in seen:
            I.pop(i)
        else:
            seen[I[i]] = 1
            i += 1
    return I

def main():
    import argparse
    from genomicode import jmath
    from genomicode import matrixlib
    from genomicode import quantnorm

    parser = argparse.ArgumentParser(
        description="Slice the rows or columns of a matrix.")
    parser.add_argument("filename", nargs="+", help="Matrices to slice.")
    parser.add_argument(
        "--read_as_csv", default=False, action="store_true",
        help="Read as a CSV file.")
    parser.add_argument(
        "--remove_comments", default=None,
        help="Remove rows that start with this character (e.g. '#')")
    parser.add_argument(
        "--clean_only", default=False, action="store_true",
        help="Only read_as_csv and remove_comments.")
    # If the user chooses an outfile, will need to implement it for
    # clean_only as well.
    #parser.add_argument(
    #    "-o", default=None, metavar="OUTFILE", dest="outfile",
    #    help="Save to this file.  By default, writes output to STDOUT.")
    parser.add_argument(
        "-l", "--log_transform", dest="log_transform", default=False,
        action="store_true",
        help="Log transform the data.")
    parser.add_argument(
        "-q", "--quantile", dest="quantile", action="store_true",
        default=None,
        help="Quantile normalize the data.")

    # --select_row_indexes   Indexes of rows to include.
    # --select_row_ids       Ids of rows to include.
    # --select_row_geneset   <gmx/gmt_file>[,<geneset>]
    # --align_row_file       <filename>
    # --add_row_annot        <gmx/gmt_file>[,<geneset>]
    # --remove_row_annot     <name>                          (multiple)
    # --select_col_indexes   Indexes of columns to include.
    # --relabel_col_ids      <gmx/gmt_file>[,<geneset>]
    # --align_col_file
    group = parser.add_argument_group(title="Column operations")
    group.add_argument(
        "--select_col_indexes", default=None,  
        help="Which columns to include e.g. 1-5,8 (1-based, inclusive)."
       )
    group.add_argument(
        "--select_col_annotation", default=None, 
        help="Include only the cols where the annotation contains a "
        "specific value.  Format: <txt_file>,<header>,<value>[,<value,...]")
    group.add_argument(
        "--select_col_genesets", default=None,
        help="Include only the samples from this geneset.  "
        "Format: <txt/gmx/gmt_file>[,<geneset>,<geneset>,...]")
    group.add_argument(
        "--remove_col_ids", default=None,  
        help="Comma-separated list of IDs to remove.")
    group.add_argument(
        "--relabel_col_ids", default=None,  
        help="Relabel the column IDs.  Format: <txt/gmx/gmt_file>,<geneset>.  "
        "One of the genesets in the file must match the current column IDs.")
    group.add_argument(
        "--align_col_file", default=None,
        help="Align the cols to this other matrix file.")
    group.add_argument(
        "--ignore_missing_cols", default=False, action="store_true",
        help="Ignore any cols that can't be found.")
    group.add_argument(
        "--filter_duplicate_cols", default=False, action="store_true",
        help="If a column is found multiple times, keep only the first one.")

    group = parser.add_argument_group(title="Row operations")
    group.add_argument(
        "--select_row_indexes", default=None,  
        help="Which rows to include e.g. 1-50,75 (1-based, inclusive)."
       )
    group.add_argument(
        "--select_row_ids", default=[], action="append",
        help="Comma-separate list of IDs to include.")
    group.add_argument(
        "--select_row_annotation", default=None, 
        help="Include only the rows where the annotation contains a "
        "specific value.  Format: <txt_file>,<header>,<value>[,<value,...]")
    group.add_argument(
        "--select_row_numeric_annotation", default=None, 
        help="Include only the rows where the annotation contains a "
        "numeric value.  Format: <txt_file>,<header>,<value>[,<value,...].  "
        'If <value> starts with a "<", then will only find the rows where '
        "the annotation is less than <value>.  "
        'The analogous constraint will be applied for ">".')
    group.add_argument(
        "--select_row_genesets", default=None,
        help="Include only the IDs from this geneset.  "
        "Format: <txt/gmx/gmt_file>[,<geneset>,<geneset>,...]")
    group.add_argument(
        "--filter_row_mean", default=None,
        help="Remove this percentage of rows that have the lowest mean.  "
        "Should be between 0 and 1.")
    group.add_argument(
        "--filter_row_var", default=None,
        help="Remove this percentage of rows that have the lowest variance.  "
        "Should be between 0 and 1.")

    group.add_argument(
        "--add_row_annot", default=None,
        help="Add a geneset as a new annotation for the matrix.  "
        "The format should be: <txt/gmx/gmt_file>,<geneset>[,<geneset>].  "
        "Each geneset in the file should contain the same number of "
        "genes as the matrix.  One of the genesets should be align-able "
        "to the IDs of this matrix.")
    group.add_argument(
        "--remove_row_annot", action="append", default=[],
        help="Remove this annotations from the matrix.")
    group.add_argument(
        "--rename_row_annot", action="append", default=[],
        help="Rename this header.  "
        "The format should be: <old_name>,<new_name>.")
    group.add_argument(
        "--move_row_annot", action="append", default=[],
        help="Move this header.  "
        "The format should be: <header>,<old_index>,<new_index>.  "
        "The indexes are 1-based.")
    group.add_argument(
        "--align_row_file", default=None,
        help="Align the rows to this other matrix file.")
    group.add_argument(
        "--ignore_missing_rows", default=False, action="store_true",
        help="Ignore any rows that can't be found.")
    
    args = parser.parse_args()
    assert len(args.filename) >= 1

    x = read_matrices(
        args.filename, args.read_as_csv, args.remove_comments, args.clean_only)
    fmt_module, matrices = x
    if len(matrices) == 1:
        MATRIX = matrices[0]
    else:
        # Merge the matrices into one big file.
        matrices = matrixlib.align_rows(*matrices)
        MATRIX = matrixlib.merge_matrices(*matrices)
    if not MATRIX.nrow():
        return

    # Slice to a submatrix.
    I1 = find_row_indexes(MATRIX, args.select_row_indexes)
    I2 = find_row_ids(MATRIX, args.select_row_ids)
    I3 = find_row_genesets(MATRIX, args.select_row_genesets)
    I4 = find_row_annotation(MATRIX, args.select_row_annotation)
    I5 = find_row_numeric_annotation(
        MATRIX, args.select_row_numeric_annotation)
    I6 = find_row_mean_var(MATRIX, args.filter_row_mean, args.filter_row_var)
    I_row = _intersect_indexes(I1, I2, I3, I4, I5, I6)
    I1 = find_col_indexes(MATRIX, args.select_col_indexes)
    I2 = find_col_genesets(MATRIX, args.select_col_genesets)
    I3 = find_col_annotation(MATRIX, args.select_col_annotation)
    I_col = _intersect_indexes(I1, I2, I3)
    MATRIX = MATRIX.matrix(I_row, I_col)

    # Remove row annotations.
    for name in args.remove_row_annot:
        MATRIX = remove_row_annot(MATRIX, name)

    # Add row annotations.
    MATRIX = add_row_annot(MATRIX, args.add_row_annot)

    # Rename the row annotations.  Do this after removing and adding.
    for x in args.rename_row_annot:
        MATRIX = rename_row_annot(MATRIX, x)

    # Move the row annotations.  Do this at the end.
    for x in args.move_row_annot:
        MATRIX = move_row_annot(MATRIX, x)

    # Relabel the column IDs.
    MATRIX = relabel_col_ids(MATRIX, args.relabel_col_ids)

    # Remove col IDs.  Do this after relabeling.
    MATRIX = remove_col_ids(MATRIX, args.remove_col_ids)

    # Filter after relabeling.
    MATRIX = remove_duplicate_cols(MATRIX, args.filter_duplicate_cols)

    # Align to the align_file.  Do this as close to the end as
    # possible, after everything else removed and added.
    I_row = align_rows(MATRIX, args.align_row_file, args.ignore_missing_rows)
    MATRIX = MATRIX.matrix(I_row, None)
    I_col = align_cols(MATRIX, args.align_col_file, args.ignore_missing_cols)
    MATRIX = MATRIX.matrix(None, I_col)

    # Log transform, if requested.
    if args.log_transform:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)

    # Quantile normalize, if requested.
    if args.quantile:
        MATRIX = quantnorm.normalize(MATRIX)

    # Write the outfile (in the same format).
    handle = sys.stdout
    #if args.outfile:
    #    handle = open(args.outfile, 'w')
    fmt_module.write(MATRIX, handle)
    

if __name__ == '__main__':
    main()
