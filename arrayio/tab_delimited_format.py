"""Generic tab-delimited format.

<ROW_HEADER_0>  [<ROW_HEADER_N>, ...]  <name_0> ...
<COL_HEADER_0>
...
<name_0>
<name_1>
...
<name_n>

If the file includes samples names, they will be stored in a special
column header SAMPLE_NAME.

Functions:
read
write
is_format
is_matrix

"""
import os, sys

SAMPLE_NAME = "_SAMPLE_NAME"

def is_format(locator_str, hrows=None, hcols=None):
    from genomicode import filelib
    if not filelib.exists(locator_str):
        return False

    # Read 5 lines and check the headers.  If the file is small, this
    # may contain fewer than 5 lines.
    handle = filelib.openfh(locator_str)
    lines = [handle.readline() for i in range(5)]
    handle.close()   # need to close it properly, or gunzip might not die.
    lines = [x for x in lines if x]

    # Make sure there's at least 1 line.
    if not lines:
        return False

    # All rows should contain a tab.
    for line in lines:
        if "\t" not in line:
            return False
    
    # All rows should contain the same number of columns.
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]
    for cols in matrix:
        if len(cols) != len(matrix[0]):
            return False
        
    return True

def is_matrix(X):
    # Any matrix can be a tab-delimited format.
    return True

## def _num_headers(matrix):
##     """Return (# row headers, # col headers)."""
##     # Try to find the number of rows and columns that contain header
##     # information.  In general, assume that numbers are data and
##     # headers contain characters.  Headers are at the beginnings of
##     # the matrix.
##     if not matrix:
##         return 0, 0
##     num_rows, num_cols = len(matrix), len(matrix[0])

##     # Make sure each row contains the same number of columns.
##     for cols in matrix:
##         assert len(cols) == len(matrix[0]), "matrix row length mismatch"

##     # First, look for rows and columns that contain characters.  Check
##     # from the end of the matrix, and move towards the front.
##     # 
##     # This algorithm is broken.  If a column is mostly numeric, but
##     # has a few non-numeric entries, then this will consider the whole
##     # column numeric because the last row is numeric.
##     last_row = matrix[num_rows-1]
##     for c in range(num_cols-1, -1, -1):
##         if not _is_numeric(last_row[c]):
##             break
##     else:
##         c = -1
##     hcols = c+1

##     last_col = [x[num_cols-1] for x in matrix]
##     for r in range(num_rows-1, -1, -1):
##         if not _is_numeric(last_col[r]):
##             break
##     else:
##         r = -1
##     hrows = r+1
##     print "HEADERS_1", hrows, hcols

##     # Check for row headers that match the pattern:
##     # <HEAD1>  <HEAD2>   <HEAD3>   <DATA>
##     # <ROW>    <blank>   <blank>  <number>
##     # Because <ROW> has <blank> value, <ROW> must be a header.  If it
##     # contains information about a specific gene, then <HEAD2> and
##     # <HEAD3> would have annotations.
##     # This will not be found by the previous filter, because it looked
##     # only in the last column, which should have values for <HEAD2>
##     # and <HEAD3>.
##     # This will fail if the row begins with a bunch of missing values.
    
##     # Exception: If the whole column is blank, and none of the genes
##     # have any annotation, then ignore this column.  This should not
##     # happen if I detect and remove blank columns above.
##     if hrows and hcols:
##         # Look for blank columns and exclude them.
##         while hrows < len(matrix) and hcols > 1:
##             all_space = 1
##             for i in range(1, hcols):
##                 if matrix[hrows][i].strip() != "":
##                     all_space = 0
##                     break
##             if not all_space:
##                 break
##             hrows += 1
##     print "HEADERS_2", hrows, hcols
            
##     # The values of some columns are numeric, but nevertheless should
##     # be considered headers.  If known headers are given, check for
##     # them.
##     if hrows:
##         col_headers = [x.upper() for x in matrix[0]]
##         while hcols < len(col_headers) and col_headers[hcols] in [
##             "GID", "NA", "ID", "NAME", "LOCUSLINK",
##             "GWEIGHT", "GORDER", "GCLUSTER"]:
##             hcols += 1
##     if hcols:
##         row_headers = [x[0].upper() for x in matrix]
##         while hrows < len(row_headers) and row_headers[hrows] in [
##             "GID", "AID", "EWEIGHT", "EORDER", "ACLUSTER"]:
##             hrows += 1
##     #print "HEADERS_3", hrows, hcols

##     # Check for row headers that match the pattern:
##     # <HEADER1>  <HEADER2>  <SAMPLE#>  <SAMPLE#>  <SAMPLE#>
##     # <NAME>     <NAME>        ...        ...        ...
##     # If the sample names are all numbers, then it might be mistaken
##     # for data.  Try to identify the case where SAMPLEs are given as
##     # integers, but the data is given as floats.
##     if hcols and not hrows:
##         samples_are_ints = False
##         for x in matrix[0][hcols:]:
##             if not _is_int(x):
##                 break
##         else:
##             samples_are_ints = True
##         data_are_float = False
##         for x in matrix[hrows+1][hcols:]:
##             if not _is_float_not_int(x):
##                 break
##         else:
##             data_are_float = True
##         if samples_are_ints and data_are_float:
##             hrows += 1
            
##     # Now check from the front of the matrix to make sure that the
##     # row and column we ended up at contains numbers.
##     while (hrows < num_rows and hcols < num_cols and
##            not _is_numeric(matrix[hrows][hcols])):
##         hcols += 1
##     #print "HEADERS_4", hrows, hcols

##     return hrows, hcols

def read(handle, hrows=None, hcols=None, datatype=float):
    import math
    from genomicode import filelib
    from genomicode import Matrix
    from genomicode import jmath
    from genomicode import iolib
    import util
    import const
    # Format:
    # - gene x experiment
    # - optional header row
    # - optional rows of sample annotations (requires header row)
    # - optional columns of gene annotations

    filename = None
    if type(handle) is type(""):
        filename = handle
    handle = filelib.openfh(handle)
    data = iolib.split_tdf(handle.read(), strip=True)

    # Sometimes people insert blank rows or columns inside the matrix.
    # Remove all of those.
    # Delete blank rows.
    i = 0
    while i < len(data):
        for x in data[i]:
            if x:
                # There's data in this row.  Ignore it.
                i += 1
                break
        else:
            del data[i]

    # Make sure each line has the same number of columns.
    num_rows = len(data)
    num_cols_all = [len(x) for x in data]
    #print num_rows, num_cols_all

    # R creates headers that contain one fewer columns than the rest
    # of the matrix, if someone requests to write out the row names.
    # Detect this case and insert a dummy header.
    all_one_fewer = True
    for i in range(1, len(data)):
        if num_cols_all[i] != num_cols_all[0]+1:
            all_one_fewer = False
            break
    # This can happen either if the length of the first row is one
    # less than every other row of the matrix, or if matrix has
    # only 1 row.  Only add a fake header if the matrix has more
    # than one row.
    if all_one_fewer and len(data) > 1:
        header_row = data[0]
        header_row.insert("ROW_NAMES", 0)
        num_cols_all[0] += 1
    
    # Make sure each line has the same number of columns.
    if num_cols_all:
        num_cols = num_cols_all[0]
        for i, nc in enumerate(num_cols_all):
            f = ""
            if filename:
                f = " [%s]" % filename
            error_msg = "Header%s has %d columns but line %d has %d." % (
                f, num_cols, i+1, nc)
            assert nc == num_cols, error_msg
    #print num_rows, num_cols; sys.exit(0)

    # Sometimes, a user might cluster a matrix with an empty column
    # using Cluster 3.0.  In this case, Cluster 3.0 will preserve the
    # empty column, except for a "1.000000" for the EWEIGHT row
    # header.  Try to detect this case and remove the "1.000000".
    last_col = [x[-1] for x in data]
    #non_empty = [x for x in last_col if x.strip()]
    non_empty = [x for x in last_col if x]
    value = None
    if len(non_empty) == 1:
        value = jmath.safe_float(non_empty[0])
    if value is not None and abs(value-1.00) < 1E-10:
        for i in range(num_rows):
            data[i][-1] = ""

    # Matlab appends blank columns to the end.  Delete columns that
    # are completely blank.
    # DEBUG.
    #handle = open("/home/jchang/debug.txt", 'w')
    #for x in data:
    #    print >>handle, "\t".join(map(str, x))
    #handle.close()
    i = 0
    while data and data[0] and i < len(data[0]):
        # Assumes that every row is the same length.
        for row in data:
            # There's data in this column.  Ignore it.
            if row[i]:
                i += 1
                break
        else:
            [x.pop(i) for x in data]

    if not data:
        return Matrix.InMemoryMatrix([])

    # If the rows and cols not explicitly specified, then try to guess
    # them from the file.
    if hrows is None or hcols is None:
        hr, hc = util.num_headers(data)
        if hrows is None:
            hrows = hr
        if hcols is None:
            hcols = hc
    #print "HEADERS", hrows, hcols
    #num_genes, num_arrays = num_rows-hrows, num_cols-hcols

    # Pull out the row names from the columns.
    row_names = {}  # header -> list of names (1 for each gene)
    row_order = []  # in-order list of the headers
    if hcols:
        if hrows:
            # If a header row is provided, then the names of these
            # annotations are provided in the header.
            row_order = data[0][:hcols]
        else:
            # No header row.  Make default name for these annotations.
            ndigits = int(math.ceil(math.log(hcols, 10)))
            row_order = ["ANNOT%*d" % (ndigits, i+1) for i in range(hcols)]
        # Strip extraneous whitespace from the header names.
        # Not necessary.  Handled now in split_tdf.
        #row_order = [x.strip() for x in row_order]
        
        # Sometimes the format detection can go wrong and a GCT file
        # will slip through to here.  If this occurs, a "duplicate
        # header" exception will be generated.  Check for this and
        # generate a more meaningful error message.
        if(row_order[0] == "#1.2" and len(row_order) > 1 and
           row_order[1] == "" and row_order[-1] == ""):
            raise AssertionError, "ERROR: It looks like a GCT file was missed."
        for i, header in enumerate(row_order):
            names = [x[i] for x in data[hrows:]]
            assert header not in row_names, "duplicate header: %s" % header
            row_names[header] = names

    # Pull out the column names.
    col_names = {}   # header -> list of names (1 for each array)
    col_order = []
    if hrows:
        for i in range(1, hrows):
            header = data[i][0]
            names = data[i][hcols:]
            assert header not in col_names, "duplicate name: %s" % header
            # Strip extraneous whitespace from the header names.
            # Not necessary.  Handled now in split_tdf.
            #header = header.strip()
            col_order.append(header)
            col_names[header] = names

    # Now extract the expression values.
    matrix = data
    if hrows or hcols:
        matrix = [x[hcols:] for x in matrix[hrows:]]

    # Pull out the sample names.
    sample_names = None
    if hrows:
        # If a header is provided, then use these as the column names.
        sample_names = data[0][hcols:]
    if sample_names:
        col_names[SAMPLE_NAME] = sample_names
        col_order.append(SAMPLE_NAME)

    if datatype is None:
        convert_fn = None   # no conversion
    elif datatype is int:
        convert_fn = jmath.safe_int
    elif datatype is float:
        convert_fn = jmath.safe_float
    else:
        # Assume that I was passed a function.
        convert_fn = datatype
        
    if convert_fn:
        check_each_row = False
        try:
            matrix = [map(convert_fn, x) for x in matrix]
        except ValueError, err1:
            if str(err1) == "empty string for float()":
                check_each_row = True
            elif str(err1).startswith("invalid literal for float()"):
                check_each_row = True
            elif str(err1).startswith("could not convert string to float"):
                check_each_row = True
            else:
                raise
        if check_each_row:
            # If there was an exception, then check each row carefully
            # to try to pinpoint the problem.
            for i, x in enumerate(matrix):
                try:
                    map(convert_fn, x)
                except ValueError, err2:
                    raise ValueError, "%s\nProblem with row %d: %s" % (
                        str(err2), i+1, data[hrows+i])
            raise AssertionError, "Error converting values."

    # Set ROW_ID and COL_ID to reasonable defaults.
    synonyms = {}
    if SAMPLE_NAME in col_names:
        synonyms[const.COL_ID] = SAMPLE_NAME
    if row_order:
        # Bug: This should be the first column with unique values.
        synonyms[const.ROW_ID] = row_order[0]

    X = Matrix.InMemoryMatrix(
        matrix, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    X = Matrix.add_synonyms(X, synonyms)
    return X

CLEAN_RE = None
CLEAN_DISALLOWED = None
def _clean(s, disallowed=None):
    # Make sure there are no disallowed characters in the string s.
    global CLEAN_RE
    global CLEAN_DISALLOWED
    
    import re

    disallowed = disallowed or "\r\n\t"
    if CLEAN_RE is None or CLEAN_DISALLOWED != disallowed:
        CLEAN_RE = re.compile("[%s]" % disallowed)
        CLEAN_DISALLOWED = disallowed
    s, count = CLEAN_RE.subn(" ", s)
    s = s.strip()
    return s

def _clean_many(l, disallowed=None):
    l = [_clean(x, disallowed=disallowed) for x in l]
    return l
    
def write(X, handle):
    from genomicode import filelib
    from genomicode import iolib

    assert is_matrix(X)
    if type(handle) is type(""):
        handle = open(handle, 'w')

    row_names = X.row_names()
    col_names = X.col_names()

    M_out = []  # Matrix to write out.

    # Print out the header row if there are row headers or sample
    # names.
    if row_names:
        header = row_names
        if SAMPLE_NAME in col_names:
            header = header + X.col_names(SAMPLE_NAME)
        M_out.append(header)
        #header = _clean_many(header)
        #print >>handle, "\t".join(header)

    # Print out the column annotations.
    for header in col_names:
        if header == SAMPLE_NAME:
            continue
        x = [header] + [""]*(len(row_names)-1) + X.col_names(header)
        M_out.append(x)
        #x = _clean_many(map(str, x))
        #print >>handle, "\t".join(x)

    # Print out the row ids and data.
    nrow, ncol = X.dim()
    M = X.slice()
    header2rownames = {}
    for header in row_names:
        header2rownames[header] = X.row_names(header)
    for i in range(nrow):
        #names = [X.row_names(header)[i] for header in row_names]
        names = [header2rownames[header][i] for header in row_names]
        # M[i] might be tuples.
        values = list(M[i])
        x = names + values
        M_out.append(x)
    iolib.cleanwrite(M_out, handle)

