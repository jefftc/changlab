"""

Functions:
align_rows                 Align multiple matrices by rows.
align_cols                 Align multiple matrices by column.
are_rows_aligned
are_cols_aligned
assert_rows_aligned
describe_unaligned_rows

align_rows_to_annot        Align a matrix to a row annotation.
align_cols_to_annot
align_rows_to_many_annots
align_cols_to_many_annots

find_row_header            Find the row header that contains specified annots.
find_col_header            Find the col header that contains specified annots.
find_best_row_header
find_best_col_header

read_matrices
merge_matrices
merge_gct_matrices

"""

# Optimization: call hash_R less

def align_rows(*matrices):
    """Aligns matrices by ROW_ID.  Return a list of the matrices after
    the rows are aligned.  Raises an exception if no rows are common
    to all matrices.

    """
    import arrayio
    import hashlib

    header = arrayio.ROW_ID
    if not matrices:
        return []
    for m in matrices:
        assert m.nrow() > 0 and m.ncol() > 0

    # Make sure each of the matrices contain the row ID.
    for m in matrices:
        assert header in m.row_names() or header in m._synonyms, \
            "Matrix is missing row IDs (%s)." % header

    # Get the intersection of the IDs.
    all_ids = [None]*len(matrices) # matrix index -> list of hids
    for i in range(len(matrices)):
        all_ids[i] = hashlib.hash_many_geneids(matrices[i].row_names(header))
    
    ids_hash = {}  # hid -> matrix_i -> row_i, id
    for i, ids in enumerate(all_ids):
        for j, id_ in enumerate(ids):
            hid = hashlib.hash_geneid(id_)
            # If I have not seen this hashed ID before, then add it to
            # the dictionary.
            if hid not in ids_hash:
                ids_hash[hid] = {}
            # If a dataset has duplicate IDs, then take the first one only.
            if i in ids_hash[hid]:
                continue
            ids_hash[hid][i] = j, id_

    # Use only the IDs that occur in all data files.  Use the order of
    # the samples from first data set.
    ids = [x for x in all_ids[0] if len(ids_hash[x]) == len(matrices)]
    assert len(ids) > 0, "The data sets all have different row IDs."
            
    # Align the rows by the ids.
    aligned = [None] * len(matrices)
    for i in range(len(matrices)):
        I = [ids_hash[x][i][0] for x in ids]
        x = matrices[i].matrix(row=I)
        aligned[i] = x

    assert are_rows_aligned(*aligned), "matrices not aligned"

    return aligned

def align_cols(*matrices):
    """Aligns matrices by COL_ID.  Return a list of the matrices after
    the columns are aligned.  Raises an exception if no columns are
    common to all matrices.

    """
    import arrayio
    import hashlib

    header = arrayio.COL_ID
    if not matrices:
        return []
    for m in matrices:
        assert m.nrow() > 0 and m.ncol() > 0

    # Get the intersection of the IDs.
    all_ids = [None]*len(matrices) # matrix index -> list of hids
    for i in range(len(matrices)):
        all_ids[i] = hashlib.hash_many_sampleids(matrices[i].col_names(header))
    
    ids_hash = {}  # hid -> matrix_i -> col_i, id
    for i, ids in enumerate(all_ids):
        for j, id_ in enumerate(ids):
            hid = hashlib.hash_sampleid(id_)
            # If I have not seen this hashed ID before, then add it to
            # the dictionary.
            if hid not in ids_hash:
                ids_hash[hid] = {}
            # If a dataset has duplicate IDs, then take the first one only.
            if i in ids_hash[hid]:
                continue
            ids_hash[hid][i] = j, id_
            
    # Use only the IDs that occur in all data files.  Use the order of
    # the samples from first data set.
    ids = [x for x in all_ids[0] if len(ids_hash[x]) == len(matrices)]
    assert len(ids) > 0, "The data sets have different column IDs."

    # Align the columns by the ids.
    aligned = [None] * len(matrices)
    for i in range(len(matrices)):
        I = [ids_hash[x][i][0] for x in ids]
        x = matrices[i].matrix(col=I)
        aligned[i] = x

    assert are_cols_aligned(*aligned), "matrices not aligned"

    return aligned

def are_rows_aligned(*matrices):
    import arrayio
    import hashlib

    header = arrayio.ROW_ID
    
    if len(matrices) <= 1:
        return True

    # Check the number of rows in each matrix.
    for i in range(1, len(matrices)):
        if matrices[0].nrow() != matrices[i].nrow():
            return False

    # Check the names of the rows.
    hnames = [hashlib.hash_many_geneids(x.row_names(header)) for x in matrices]
    for i in range(1, len(matrices)):
        if hnames[0] != hnames[i]:
            return False
    return True

def are_cols_aligned(*matrices):
    import arrayio
    import hashlib

    header = arrayio.COL_ID

    if len(matrices) <= 1:
        return True

    # Check the number of columns in each matrix.
    for i in range(1, len(matrices)):
        if matrices[0].ncol() != matrices[i].ncol():
            return False

    # Check the names of the columns.
    hnames = [
        hashlib.hash_many_sampleids(x.col_names(header)) for x in matrices]
    for i in range(1, len(matrices)):
        if hnames[0] != hnames[i]:
            return False
    return True

def assert_rows_aligned(*matrices):
    # Ignore None.
    matrices = [x for x in matrices if x]
    assert are_rows_aligned(*matrices)

def describe_unaligned_rows(*matrices):
    # Return a text string that describes where the rows are not aligned.
    import arrayio
    import parselib
    import hashlib

    header = arrayio.ROW_ID
    
    if len(matrices) <= 1:
        return "Only 1 matrix.  Must be aligned."

    # Check the number of rows in each matrix.
    num_rows = [x.nrow() for x in matrices]
    if min(num_rows) != max(num_rows):
        x = "Matrices have differing number of rows: %s" % ", ".join(
            map(parselib.pretty_int, num_rows))
        return x

    # Check the names of the rows.
    names = [x.row_names(header) for x in matrices]
    hnames = [hashlib.hash_many_geneids(x.row_names(header)) for x in matrices]
    bad_rows = []
    for i in range(matrices[0].nrow()):
        unaligned =False
        for j in range(1, len(matrices)):
            if hnames[0][i] != hnames[j][i]:
                unaligned = True
        if unaligned:
            x = [names[j][i] for j in range(len(matrices))]
            x = "Row %s: %s" % (parselib.pretty_int(i+1), ", ".join(x))
            bad_rows.append(x)
    if not bad_rows:
        return "Matrices are aligned."

    total_bad = len(bad_rows)
    if total_bad > 10:
        bad_rows = bad_rows[:10]
        bad_rows.append("...")
    x = "%s of %s rows are unaligned." % (
        parselib.pretty_int(total_bad),
        parselib.pretty_int(matrices[0].nrow()))
    lines = [x] + bad_rows
    return "\n".join(lines)

def align_rows_to_annot(
        MATRIX, row_names, header=None, hash=False, get_indexes=False):
    # Return tuple of aligned (MATRIX, row_names).  If get_indexes is
    # True, then will return a tuple (I_MATRIX, I_row_names) of the
    # indexes to align the matrices instead.  If there are no common
    # names, then MATRIX and row_names will be empty.  If header is
    # given, then will use the row names from that header.  Otherwise,
    # will search through all possible headers for the best match.
    # XXX document hash
    x = _align_to_annot_I(MATRIX, row_names, header, hash, True)
    I1, I2 = x
    if get_indexes:
        return I1, I2
    MATRIX = MATRIX.matrix(I1, None)
    row_names = [row_names[x] for x in I2]
    return MATRIX, row_names

def align_cols_to_annot(
        MATRIX, col_names, header=None, hash=False, get_indexes=False):
    # Return tuple of aligned (MATRIX, col_names).
    x = _align_to_annot_I(MATRIX, col_names, header, hash, False)
    I1, I2 = x
    if get_indexes:
        return I1, I2
    MATRIX = MATRIX.matrix(None, I1)
    col_names = [col_names[x] for x in I2]
    return MATRIX, col_names

def align_rows_to_many_annots(
        MATRIX, many_row_names, header=None, hash=False, get_indexes=False):
    # Return tuple of aligned (MATRIX, row_names).  If get_indexes is
    # True, then return (I_MATRIX, I_row_names, index into
    # many_row_names).
    x = _align_to_many_annots_I(MATRIX, many_row_names, header, hash, True)
    I_matrix, I_names, index = x
    if get_indexes:
        return I_matrix, I_names, index
    
    MATRIX = MATRIX.matrix(I_matrix, None)
    row_names = many_row_names[index]
    row_names = [row_names[x] for x in I_names]
    return MATRIX, row_names

def align_cols_to_many_annots(
        MATRIX, many_col_names, header=None, hash=False, get_indexes=False):
    # many_col_names is list of col_names, where col_names is a list
    # of the column names to be matched to the MATRIX.  Return tuple
    # of aligned (MATRIX, col_names).  If get_indexes is True, then
    # return (I_MATRIX, I_col_names, index into many_col_names).  XXX
    # BEST MATCH?
    x = _align_to_many_annots_I(MATRIX, many_col_names, header, hash, False)
    I_matrix, I_names, index = x
    if get_indexes:
        return I_matrix, I_names, index
    
    MATRIX = MATRIX.matrix(I_matrix, None)
    col_names = many_col_names[index]
    col_names = [col_names[x] for x in I_names]
    return MATRIX, col_names

def _align_to_many_annots_I(MATRIX, many_names, header, hash, is_row):
    best_I_matrix = best_I_names = best_index = None
    for i, names in enumerate(many_names):
        x = _align_to_annot_I(MATRIX, names, header, hash, is_row)
        I_matrix, I_names = x
        #print names, header, I_matrix
        if best_I_matrix is None or len(I_matrix) > len(best_I_matrix):
            best_I_matrix = I_matrix
            best_I_names = I_names
            best_index = i
    return best_I_matrix, best_I_names, best_index

def _align_to_annot_I(MATRIX, names, header, hash, is_row):
    get_names = MATRIX.row_names
    if not is_row:
        get_names = MATRIX.col_names
        
    headers = [header]
    if header is None:
        headers = get_names()
    assert headers
    
    best_I_matrix = best_I_names = None
    for header in headers:
        x = _align_header_to_annot_I(MATRIX, names, header, hash, is_row)
        I_matrix, I_names = x
        if best_I_matrix is None or len(I_matrix) > len(best_I_matrix):
            best_I_matrix, best_I_names = I_matrix, I_names
    return best_I_matrix, best_I_names

def _align_header_to_annot_I(MATRIX, names, header, hash, is_row):
    import jmath
    import hashlib

    get_names = MATRIX.row_names
    if not is_row:
        get_names = MATRIX.col_names
    
    matrix_names = get_names(header)

    h_names = names
    h_matrix_names = matrix_names
    if hash:
        #h_names = [hashlib.hash_R(x) for x in names]
        #h_matrix_names = [hashlib.hash_R(x) for x in matrix_names]
        h_names = hashlib.hash_R_many(names)
        h_matrix_names = hashlib.hash_R_many(matrix_names)
    
    I_names = jmath.match(h_matrix_names, h_names)
    I_MATRIX = []
    for i in range(len(I_names)):
        if I_names[i] is not None:
            I_MATRIX.append(i)
    I_names = [x for x in I_names if x is not None]
    assert len(I_MATRIX) == len(I_names)
    return I_MATRIX, I_names
    
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

def find_row_header(MATRIX, row_names, hash=False):
    # Find the column (row header) from the MATRIX that contains all
    # the row_names.  If multiple headers match, then just return the
    # first one.  If none match, return None.
    x = find_best_row_header(MATRIX, row_names, hash=hash)
    header, num_matches, found, missing = x
    if num_matches == len(row_names):
        return header
    return None

def find_col_header(MATRIX, col_names, hash=False):
    x = best_col_header(MATRIX, col_names, hash=hash)
    header, num_matches, found, missing = x
    if num_matches == len(col_names):
        return header
    return None

def _find_best_header(MATRIX, names, hash, is_row):
    import hashlib

    get_names = MATRIX.row_names
    if not is_row:
        get_names = MATRIX.col_names

    h_names = names
    if hash:
        # Hash the names for comparison.
        #h_names = [hashlib.hash_R(x) for x in names]
        h_names = hashlib.hash_R_many(names)

    # Count the number of names that are found in each header.
    header2found = {}    # header -> list of names found
    header2missing = {}  # header -> list of names not found
    for header in get_names():
        matrix_names = get_names(header)
        h_matrix_names = matrix_names
        if hash:
            h_matrix_names = [hashlib.hash_R(x) for x in matrix_names]
        h_matrix_names = {}.fromkeys(h_matrix_names)

        x1 = [x for x in h_names if x in h_matrix_names]
        x2 = [x for x in h_names if x not in h_matrix_names]
        header2found[header] = x1
        header2missing[header] = x2

    # Find the best match.
    best_header = best_found = best_num_missing = missing = None
    for header in get_names():
        found = header2found[header]
        missing = header2missing[header]
        num_missing = len(missing)
        if best_header is None or num_missing < best_num_missing:
            best_header = header
            best_found = found
            best_num_missing = num_missing
            best_missing = missing
    return best_header, len(names)-best_num_missing, best_found, best_missing

def find_best_row_header(MATRIX, row_names, hash=False):
    # Find the column (row header) from the MATRIX that contains the
    # most row_names.  If multiple headers match the same number, then
    # just return the first one.  Return a tuple of (header,
    # num_matches, list of matches, list of mismatches).
    return _find_best_header(MATRIX, row_names, hash, True)

def find_best_col_header(MATRIX, col_names, hash=False):
    return _find_best_header(MATRIX, col_names, hash, False)

def read_matrices(filenames, cache=None):
    """Read a list of matrices and align them.  filenames is a list of
    the matrix files to read.  Returns a tuple where the first element
    is a list of the matrices read, and the second is the aligned
    matrix.

    cache is an optional dictionary of filename to matrix.  This can
    be used to prevent re-reading of matrices.

    """
    import copy
    import arrayio
    import filelib

    for filename in filenames:
        assert filelib.exists(filename), "File not found: %s" % filename

    # Load the files.
    DATA = []
    for filename in filenames:
        if cache is not None and filename in cache:
            x = copy.deepcopy(cache[filename])
        else:
            try:
                x = arrayio.read(filename)
            except (SystemError, KeyboardInterrupt, MemoryError), x:
                raise
            except Exception, x:
                # Can diagnose which file failed here.
                # raise
                raise Exception, "Problem reading %s: %s" % (
                    repr(filename), str(x))
            if cache is not None:
                cache[filename] = x
        DATA.append(x)

    #for d, filename in zip(DATA, filenames):
    #    f = os.path.split(filename)[1]
    #    print "%s has %d genes and %d samples." % (f, d.nrow(), d.ncol())

    # Align the matrices.
    ALIGNED = align_rows(*DATA)
    #if DATA:
    #    print "Merged file has %d genes." % DATA[0].nrow()
    #sys.stdout.flush()

    return DATA, ALIGNED

def merge_matrices(*matrices):
    import arrayio
    import Matrix

    assert len(matrices)
    assert are_rows_aligned(*matrices)

    # Get the values in list of lists format.
    X = [x.value() for x in matrices]
    X_all = []
    for i in range(len(X[0])):   # each row
        x = []
        for j in range(len(X)):  # each matrix
            x.extend(X[j][i])
        X_all.append(x)

    # The sample names may not be unique.
    matrix0 = matrices[0]
    row_order = []
    col_order = []
    row_names = {}
    col_names = {}
    synonyms = {}

    # Use the annotations from the first matrix.
    for name in matrix0.row_names():
        row_names[name] = matrix0.row_names(name)
    row_order = matrix0._row_order
    # Add any annotations that do not occur in the first matrix.
    for m in matrices[1:]:
        for name in m.row_names():
            if name in row_names:
                continue
            row_names[name] = m.row_names(name)
            row_order.append(name)

    # NotImplemented: copy over the column annotations.

    # Copy over the sample names.
    col_order = [matrix0._col_order[0]]
    x = []
    for m in matrices:
        x.extend(m.col_names(arrayio.COL_ID))
    col_names[col_order[0]] = x
    
    synonyms[arrayio.ROW_ID] = matrix0._synonyms[arrayio.ROW_ID]
    synonyms[arrayio.COL_ID] = col_order[0]

    DATA = Matrix.InMemoryMatrix(
        X_all, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #DATA = Matrix.add_synonyms(DATA, synonyms)
    return DATA

def merge_gct_matrices(*matrices):
    import arrayio
    import Matrix

    assert len(matrices)
    assert are_rows_aligned(*matrices)
    for m in matrices:
        assert arrayio.gct_format.is_matrix(m)

    # Get the values in list of lists format.
    X = [x.value() for x in matrices]
    X_all = []
    for i in range(len(X[0])):   # each row
        x = []
        for j in range(len(X)):  # each matrix
            x.extend(X[j][i])
        X_all.append(x)

    # The sample names may not be unique.
    row_order = ["NAME", "DESCRIPTION"]
    col_order = []
    row_names = {}
    col_names = {}
    synonyms = {}

    # Should be in GCT format.
    matrix0 = matrices[0]
    assert len(matrix0.row_names()) == 2, "Invalid file format."
    # Just assume the first column is the NAME and second is the
    # DESCRIPTION.  Allow more flexibility in the actual column
    # headers.
    name, description = matrix0.row_names()
    row_names["NAME"] = matrix0.row_names(name)
    row_names["DESCRIPTION"] = matrix0.row_names(description)

    for m in matrices:
        assert len(m.col_names()) == 1
    col_order = matrix0.col_names()
    assert len(col_order) == 1

    x = []
    for m in matrices:
        x.extend(m.col_names(m.col_names()[0]))
    col_names[col_order[0]] = x

    synonyms[arrayio.ROW_ID] = "NAME"
    synonyms[arrayio.COL_ID] = col_order[0]

    DATA = Matrix.InMemoryMatrix(
        X_all, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #DATA = Matrix.add_synonyms(DATA, synonyms)
    assert arrayio.gct_format.is_matrix(DATA)
    return DATA

