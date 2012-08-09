"""Implement GCT Gene Cluster Text array format.

#1.2
<# rows>  <# samples>
NAME      DESCRIPTION     <Sample>    ...
...

* NAME and DESCRIPTION are case insensitive.
* Relax this constraint and allow NAME and DESCRIPTION to be other things.

Samples names are stored in a special column header
tab_delimited_format.SAMPLE_NAME.

Functions:
read
write
is_format
is_matrix

"""
import os

def is_format(locator_str, hrows=None, hcols=None):
    from genomicode import filelib
    if not filelib.exists(locator_str):
        return False

    if hrows not in [None, 1]:
        return False
    if hcols not in [None, 2]:
        return False

    # Read 5 lines and check the headers.
    handle = filelib.openfh(locator_str)
    lines = [handle.readline() for i in range(5)]
    handle.close()   # need to close it properly, or gunzip might not die.

    lines = [x for x in lines if x]
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]

    if len(matrix) < 3:
        return False
    # First line could be just one column, or could be many columns.
    if len(matrix[0]) < 1:
        return False
    # Second line must have at least 2 columns.
    if len(matrix[1]) < 2:
        return False
    if matrix[0][0] != "#1.2":
        return False
    #if matrix[2][0].strip().upper() != "NAME":
    #    return False
    #if matrix[2][1].strip().upper() != "DESCRIPTION":
    #    return False
    return True

DIAGNOSIS = ""
def is_matrix(X):
    global DIAGNOSIS
    import tab_delimited_format as tdf

    DIAGNOSIS = ""

    if not hasattr(X, "col_names") or not hasattr(X, "row_names"):
        DIAGNOSIS = "No annotations."
        return False
    
    if tdf.SAMPLE_NAME not in X.col_names():
        DIAGNOSIS = "No samples."
        return False
    if len(X.col_names()) != 1:
        DIAGNOSIS = "Extra sample annotations."
        return False
    if len(X.row_names()) != 2:
        DIAGNOSIS = "Row annotations not right."
        return False
    # Make sure "NAME" and "DESCRIPTION" are present somewhere.
    x = [x.upper() for x in X.row_names() + X._synonyms.keys()]
    if "NAME" not in x or "DESCRIPTION" not in x:
        DIAGNOSIS = "Missing NAME and/or DESCRIPTION headers."
        return False
    #x = [x.upper() for x in X.row_headers()]
    #if sorted(x) != sorted(["NAME", "DESCRIPTION"]):
    #    return False
    return True

def read(handle, hrows=None, hcols=None, datatype=float):
    from genomicode import Matrix
    import const
    import tab_delimited_format
    from genomicode import filelib

    assert hrows is None or hrows == 1
    assert hcols is None or hcols == 2

    handle = filelib.openfh(handle)
    assert handle.readline().strip() == "#1.2"
    x = handle.readline().rstrip("\r\n").split("\t")
    assert len(x) >= 2
    num_genes, num_samples = map(int, x[:2])

    X = tab_delimited_format.read(handle, hrows=1, hcols=2, datatype=datatype)
    assert X.dim() == (num_genes, num_samples), (
        "Matrix size mismatch.\n"
        "The GCT headers indicate a matrix with %d rows and %d columns.\n"
        "However, I found %d rows and %d columns." % (
            num_genes, num_samples, X.nrow(), X.ncol()))

    #assert X.row_headers()[0].upper() == "NAME"
    #assert X.row_headers()[1].upper() == "DESCRIPTION"
    header0, header1 = X.row_names()[:2]
    synonyms = {}
    NAME, DESCRIPTION = "NAME", "DESCRIPTION"
    if header0 != NAME:
        synonyms[NAME] = header0
    if header1 != DESCRIPTION:
        synonyms[DESCRIPTION] = header1
    synonyms[const.ROW_ID] = header0
    X._synonyms.update(synonyms)
    #X = Matrix.add_synonyms(X, synonyms)
    assert is_matrix(X)

    # The GCT File Format description at the Broad Institute does not
    # require the NAMEs to be unique.
    ## Make sure the NAMEs are unique.
    #seen = {}
    #dups = {}
    #for name in X.row_annots(NAME):
    #    if name in seen:
    #        dups[name] = 1
    #    seen[name] = 1
    #dups = sorted(dups)
    #assert len(dups) < 5, "%s column has %d duplicated names." % (
    #    header0, len(dups))
    #assert not dups, "%s column has duplicated names: %s" % (
    #    header0, dups)
    
    return X

def write(X, handle):
    import tab_delimited_format

    assert is_matrix(X)
    num_genes, num_samples = X.dim()
    handle.write("#1.2\n")
    handle.write("%d\t%d\n" % (num_genes, num_samples))
    tab_delimited_format.write(X, handle)
