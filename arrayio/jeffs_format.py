"""Jeff's standard format for microarray data.

Probe.Set.ID  Description  LocusLink  Gene.Symbol  [<SAMPLE> ...]

Samples names are stored in a special column header
tab_delimited_format.SAMPLE_NAME.


Functions:
read
write
is_format
is_matrix

"""
import os
import const

ROW_HEADERS = ["Probe.Set.ID", "Description", "LocusLink", "Gene.Symbol"]
MYNAME_TO_STDNAME = [
    ("Probe.Set.ID", const.ROW_ID),
    ("Probe.Set.ID", const.AFFY_PROBESET_ID),
    ("Description", const.GENE_DESCRIPTION),
    ("LocusLink", const.GENE_ID),
    ("Gene.Symbol", const.GENE_SYMBOL),
    ]

def is_format(locator_str, hrows=None, hcols=None):
    from genomicode import filelib
    import util
    
    if hrows not in [None, 1]:
        return False
    if hcols not in [None, 4]:
        return False

    if not filelib.exists(locator_str):
        # This will only work if locator_str is a string.
        return False
    
    # Read 5 lines and check the headers.  If the file is small, this
    # may contain fewer than 5 lines.
    handle = filelib.openfh(locator_str)
    lines = [handle.readline() for i in range(5)]
    handle.close()   # need to close it properly, or gunzip might not die.
    lines = [x for x in lines if x]
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]

    # Make sure there's at least 1 line.
    if not matrix:
        return False

    header = matrix[0]
    if header[:len(ROW_HEADERS)] != ROW_HEADERS:
        return False

    # Check if there's extraneous stuff.
    nr, nc = util.num_headers(matrix)
    if nc > 4:
        return False
    
    return True
    
    #handle = filelib.openfh(locator_str)
    #x = handle.readline()
    #handle.close()   # need to close it properly, or gunzip might not die.
    #row = x.rstrip("\r\n").split("\t")
    #if row[:len(ROW_HEADERS)] == ROW_HEADERS:
    #    return True
    #return False

def is_matrix(X):
    import tab_delimited_format as tdf
    
    if not hasattr(X, "row_names") or not hasattr(X, "col_names"):
        return False
    # Should only include SAMPLE_NAME.
    if len(X.col_names()) != 1:
        return False
    if len(X.row_names()) != len(ROW_HEADERS):
        return False
    for header in X.row_names():
        if header not in ROW_HEADERS:
            return False
    return True

def read(handle, hrows=None, hcols=None, datatype=float):
    from genomicode import Matrix
    import tab_delimited_format as tdf

    assert hrows is None or hrows == 1
    assert hcols is None or hcols == 4
    
    X = tdf.read(handle, hrows=1, hcols=4, datatype=datatype)
    synonyms = {}
    for myname, stdname in MYNAME_TO_STDNAME:
        synonyms[stdname] = myname
    #X = Matrix.add_synonyms(X, synonyms)
    X._synonyms.update(synonyms)
    assert is_matrix(X)
    return X

def write(X, handle):
    import tab_delimited_format
    
    assert is_matrix(X)
    tab_delimited_format.write(X, handle)
