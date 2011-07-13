"""Implement Mike Eisen's data format for storing microarray data in
Cluster program.

<ID>   NAME  GWEIGHT  GORDER  [...]
EWEIGHT
EORDER
[...]

If the file includes samples names, they will be stored in a special
column header tab_delimited_format.SAMPLE_NAME.

Functions:
read
write
is_format
is_matrix

"""
import os

# GORDER/EORDER seems to be optional.
ROW_HEADERS = ["NAME", "GWEIGHT", "GORDER"]
COL_HEADERS = ["EWEIGHT", "EORDER"]

def is_format(locator_str):
    from genomicode import filelib
    import tab_delimited_format
    if not filelib.exists(locator_str):
        return False

    # Read 5 lines and count the headers.
    handle = filelib.openfh(locator_str)
    lines = [handle.readline() for i in range(5)]
    handle.close()   # need to close it properly, or gunzip might not die.
    lines = [x for x in lines if x]
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]

    # Make sure there's at least 1 line.
    if not matrix:
        return False

    # Has to have at least a header.
    if len(matrix) < 1:
        return False
    # All rows should contain the same number of columns.
    for cols in matrix:
        if len(cols) != len(matrix[0]):
            return False
        
    nrow, ncol = tab_delimited_format._num_headers(matrix)
    # tab_delimited_format sometimes mistakenly believes the number of
    # rows is 0 if the sample names are numbers.
    if nrow == 0 and ncol == 0:
        return False
    nrow = max(nrow, 1)
    if nrow < 1 or nrow > 3:
        return False
    # PCL format has at most 2 columns.
    if ncol > 2:
        return False
    #if ncol > 2:
    #    ncol = 2
    #if ncol < 2 or ncol > 4:
    #    return False
    assert len(matrix) >= 1
    header_def = [
        (0, 1, "NAME"),   (0, 2, "GWEIGHT"),   (0, 3, "GORDER"),
        (1, 0, "EWEIGHT"),
        (2, 0, "EORDER"),
        ]
    for row, col, name in header_def:
        if nrow > row and ncol > col:
            if matrix[row][col].strip().upper() != name:
                return False
    return True

def is_matrix(X):
    import tab_delimited_format as tdf

    if not hasattr(X, "col_names") or not hasattr(X, "row_names"):
        return False
    # Needs at least the ID header.
    if len(X.row_names()) < 1:
        return False
    if tdf.SAMPLE_NAME not in X.col_names():
        return False

    for i, name in enumerate(X.row_names()):
        if i == 0:
            # Ignore the <ID> column.  Can be named anything.
            continue
        if name not in ROW_HEADERS:
            return False
    for name in X.col_names():
        # Ignore the sample name header.
        if name == tdf.SAMPLE_NAME:
            continue
        if name not in COL_HEADERS:
            return False
    return True

def read(handle, datatype=float):
    import StringIO
    import tab_delimited_format
    from genomicode import filelib

    # Figure out the number of headers for tab_delimited_format.  If
    # sample names are numbers, then tab_delimited_format might
    # mistake the first row(s) for non-headers.
    s = filelib.openfh(handle).read()

    # Read 5 lines and check the headers.  If the file is small, this
    # may contain fewer than 5 lines.
    handle = StringIO.StringIO(s)
    lines = [handle.readline() for i in range(5)]
    lines = [x for x in lines if x]
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]

    assert len(matrix) >= 1
    assert len(matrix[0]) >= 2

    hrows, hcols = 1, 1
    if len(matrix[0]) >= 2 and matrix[0][1].strip().upper() in ROW_HEADERS:
        hcols += 1
    if len(matrix[0]) >= 3 and matrix[0][2].strip().upper() in ROW_HEADERS:
        hcols += 1
    if len(matrix[0]) >= 4 and matrix[0][3].strip().upper() in ROW_HEADERS:
        hcols += 1
    if len(matrix) >= 2 and matrix[1][0].strip().upper() in COL_HEADERS:
        hrows += 1
    if len(matrix) >= 3 and matrix[2][0].strip().upper() in COL_HEADERS:
        hrows += 1

    handle = StringIO.StringIO(s)
    X = tab_delimited_format.read(
        handle, hrows=hrows, hcols=hcols, datatype=datatype)
    assert is_matrix(X)
    return X

def write(X, handle):
    import tab_delimited_format
    
    assert is_matrix(X)
    tab_delimited_format.write(X, handle)
