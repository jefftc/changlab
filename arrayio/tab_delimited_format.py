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
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]
    matrix = _clean_tdf(matrix)

    # Make sure there's at least 1 line.
    if not matrix:
        return False

    # All rows should contain at least one column.
    for x in matrix:
        if not x:
            return False

    # All rows should contain the same number of columns.
    for x in matrix:
        if len(x) != len(matrix[0]):
            return False

    return True


def is_matrix(X):
    # Any matrix can be a tab-delimited format.
    return True


def _clean_tdf(matrix):
    # Return a cleaned up matrix.
    from genomicode import jmath
    
    # Sometimes people insert blank rows or columns inside the matrix.
    # Remove all of those.
    # Delete blank rows.
    matrix = [x for x in matrix if x]

    # The first line can contain one fewer column than the rest of the
    # matrix for two reasons.
    # 1.  R writes col names that contain one fewer column than the
    #     rest of the matrix, if row names are also requested.
    # 2.  For GEO data sets, each non-header row ends with a '\t', so
    #     they have an extra blank column.
    # Detect if either of these are the case.  If Case #1, then add a
    # dummy header.  If Case #2, then delete the extra blank columns.
    
    all_one_fewer = True
    for i in range(1, len(matrix)):
        if len(matrix[i]) != len(matrix[0])+1:
            all_one_fewer = False
            break
    last_column_blank = True
    for i in range(1, len(matrix)):
        if matrix[i][-1] != "":
            last_column_blank = False
            break
    # This can happen either if the length of the first row is one
    # less than every other row of the matrix, or if matrix has
    # only 1 row.
    if all_one_fewer and len(matrix) > 1:
        if last_column_blank:
            for i in range(1, len(matrix)):
                matrix[i] = matrix[i][:-1]
        else:
            matrix[0].insert("ROW_NAMES", 0)

    ## # Make sure each line has the same number of columns.
    ## for i, x in enumerate(matrix):
    ##     f = ""
    ##     #if filename:
    ##     #    f = " [%s]" % filename
    ##     error_msg = "Header%s has %d columns but line %d has %d." % (
    ##         f, num_cols, i + 1, nc)
    ##     assert nc == num_cols, error_msg

    # Matlab appends blank columns to the end.  Delete columns that
    # are completely blank.
    # DEBUG.
    #handle = open("/home/jchang/debug.txt", 'w')
    #for x in data:
    #    print >>handle, "\t".join(map(str, x))
    #handle.close()
    i = 0
    while matrix and matrix[0] and i < len(matrix[0]):
        for row in matrix:
            # This row is too short, or there's data in this column.
            # Ignore it.
            if i >= len(row) or row[i]:
                i += 1
                break
        else:
            [x.pop(i) for x in matrix]

    # Sometimes, a user might cluster a matrix with an empty column
    # using Cluster 3.0.  In this case, Cluster 3.0 will preserve the
    # empty column, except for a "1.000000" for the EWEIGHT row
    # header.  Try to detect this case and remove the "1.000000".
    last_col = [x[-1] for x in matrix]
    #non_empty = [x for x in last_col if x.strip()]
    non_empty = [x for x in last_col if x]
    value = None
    if len(non_empty) == 1:
        value = non_empty[0]
        try:
            value = float(value)
        except ValueError, x:
            pass
    if value is not None and \
           type(value) is type(0.0) and abs(value - 1.00) < 1E-10:
        for i in range(len(matrix)):
            matrix[i][-1] = ""

    # Strip whitespace.
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = matrix[i][j].strip()

    return matrix


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
    data = filelib.read_all_cols(handle)
    #data = [x for x in filelib.read_cols(handle)]
    #x = handle.read()
    #data = iolib.split_tdf(x, strip=True)
    #handle = filelib.read_cols(handle)
    #data = [handle.next() for i in range(100)]
    data = _clean_tdf(data)

    num_cols = len(data[0])
    for i, x in enumerate(data):
        nc = len(data[i])
        f = ""
        if filename:
            f = " [%s]" % filename
        error_msg = "Header%s has %d columns but line %d has %d." % (
            f, num_cols, i + 1, nc)
        assert nc == num_cols, error_msg
    if not data:
        return Matrix.InMemoryMatrix([])

    # If the rows and cols not explicitly specified, then try to guess
    # them from the file.
    #print "HEADERS 1", hrows, hcols
    if hrows is None or hcols is None:
        hr, hc = util.num_headers(data)
        if hrows is None:
            hrows = hr
        if hcols is None:
            hcols = hc
    #print "HEADERS 2", hrows, hcols
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
            row_order = ["ANNOT%*d" % (ndigits, i + 1) for i in range(hcols)]
        # Strip extraneous whitespace from the header names.
        # Not necessary.  Handled now in split_tdf.
        #row_order = [x.strip() for x in row_order]

        # Sometimes the format detection can go wrong and a GCT file
        # will slip through to here.  If this occurs, a "duplicate
        # header" exception will be generated.  Check for this and
        # generate a more meaningful error message.
        if(row_order[0] == "#1.2" and len(row_order) > 1 and
           row_order[1] == "" and row_order[-1] == ""):
            raise AssertionError("ERROR: It looks like a GCT file was missed.")
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
        col_order.insert(0, SAMPLE_NAME)

    if datatype is None:
        convert_fn = None   # no conversion
    elif datatype is int:
        convert_fn = jmath.safe_int
    elif datatype is float:
        convert_fn = jmath.safe_float
    else:
        # Assume that I was passed a function.
        convert_fn = datatype

    if convert_fn == jmath.safe_float:
        # Try and convert to an integer instead.
        is_int = True
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if not jmath.is_int(matrix[i][j]):
                    is_int = False
                    break
            if not is_int:
                break
        if is_int:
            convert_fn = jmath.safe_int

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
                    row = data[hrows + i]
                    raise ValueError("%s\nProblem with row %d: %s" % (
                        str(err2), i + 1, row))
            raise AssertionError("Error converting values.")

        
        


    # Set ROW_ID and COL_ID to reasonable defaults.
    synonyms = {}
    if SAMPLE_NAME in col_names:
        synonyms[const.COL_ID] = SAMPLE_NAME
    if row_order:
        # Bug: This should be the first column with unique values.
        synonyms[const.ROW_ID] = row_order[0]

    X = Matrix.InMemoryMatrix(
        matrix, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #X = Matrix.add_synonyms(X, synonyms)
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
    from genomicode import iolib

    assert is_matrix(X)
    if type(handle) is type(""):
        handle = open(handle, 'w')

    row_names = X.row_names()
    col_names = X.col_names()

    M_out = []  # Matrix to write out.

    # Print out the header row if there are row headers or sample
    # names.
    header = []
    if SAMPLE_NAME in col_names:
        header = X.col_names(SAMPLE_NAME)
    if row_names:
        header = row_names + header
    if header:
        M_out.append(header)
    #if row_names:
    #    header = row_names
    #    if SAMPLE_NAME in col_names:
    #        header = header + X.col_names(SAMPLE_NAME)
    #    M_out.append(header)
    #    #header = _clean_many(header)
    #    #print >>handle, "\t".join(header)

    # Print out the column annotations.
    for header in col_names:
        if header == SAMPLE_NAME:
            continue

        # If there are no row_names, then there is no room for column
        # names.  Skip it.
        x = X.col_names(header)
        if row_names:
            x = [header] + [""] * (len(row_names) - 1) + x
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
