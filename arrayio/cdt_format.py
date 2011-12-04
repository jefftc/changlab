"""Implement Mike Eisen's data format for storing microarray data in
Cluster program.

GID    <ID>   NAME  GWEIGHT  GORDER  [...]
AID
EWEIGHT
EORDER
[...]


Samples names are stored in a special column header
tab_delimited_format.SAMPLE_NAME.


Functions:
read
write
is_format
is_matrix

"""
import os

# GORDER/EORDER seems to be optional.
ROW_HEADERS = ["GID", "NAME", "GWEIGHT", "GORDER"]
COL_HEADERS = ["AID", "EWEIGHT", "EORDER"]

def is_format(locator_str, hrows=None, hcols=None):
    from genomicode import filelib
    import util
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
    
    # All rows should contain the same number of columns.
    for cols in matrix:
        if len(cols) != len(matrix[0]):
            return False

    nr, nc = util.num_headers(matrix)
    nrow = hrows or nr
    ncol = hcols or nc
    
    if nrow < 1 or nrow > 4:
        return False
    if ncol < 1 or ncol > 5:
        return False
    header_def = [
        (0, 0, "GID"), (0, 2, "NAME"),   (0, 3, "GWEIGHT"),   (0, 4, "GORDER"),
        (1, 0, "AID"),
        (2, 0, "EWEIGHT"),
        (3, 0, "EORDER"),
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
    if tdf.SAMPLE_NAME not in X.col_names():
        return False

    for i, name in enumerate(X.row_names()):
        if i == 1:
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

def read(handle, nrows=None, hcols=None, datatype=float):
    from genomicode import Matrix
    import const
    import tab_delimited_format

    X = tab_delimited_format.read(
        handle, hrows=hrows, hcols=hcols, datatype=datatype)
    # Set const.ROW_ID to be the <ID> column.
    assert len(X.row_names()) >= 2
    synonyms = {}
    synonyms[const.ROW_ID] = X.row_names()[1]
    X = Matrix.add_synonyms(X, synonyms)
    assert is_matrix(X)
    return X

def write(X, handle):
    import tab_delimited_format
    
    assert is_matrix(X)
    tab_delimited_format.write(X, handle)
