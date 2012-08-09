"""Implement RES format from GenePattern.

Format described at:
http://www.broad.mit.edu/cancer/software/genepattern/tutorial/gp_fileformats.html

Description    Accession    <Sample>  <EMPTY>    [...]
<EMPTY>        <desc1>                <desc2>    [...]
<num rows>
<description>  <accession>  <signal>  [AP]       [...]
[...]

* Each Accession must be unique.
* The columns for Description and Accession can be switched.
* Format of the sample descriptions are:
  <description>[/scale factor=<scale>]
* The second line contains one fewer columns than the first line.
* The blank column contains absent or present calls [AP].
* The GenePattern PreprocessDataset module generates res format files
  where the last column is empty.



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

    # Read 5 lines and count the headers.
    handle = filelib.openfh(locator_str)
    lines = [handle.readline() for i in range(5)]
    handle.close()   # need to close it properly, or gunzip might not die.
    lines = [x for x in lines if x]
    matrix = [line.rstrip("\r\n").split("\t") for line in lines]

    if len(matrix) < 3:
        return False

    # Line 3 should contain only 1 column.
    if len(matrix[2]) != 1:
        return False

    # Line 1 contains 1 more column than line 2.
    if len(matrix[0]) != len(matrix[1])+1:
        return False

    if len(matrix[0]) < 2:
        return False
    x = [x.upper() for x in matrix[0][:2]]
    if sorted(x) != sorted(["ACCESSION", "DESCRIPTION"]):
        return False

    return True

DIAGNOSIS = ""
def is_matrix(X):
    global DIAGNOSIS
    import tab_delimited_format as tdf
    
    DIAGNOSIS = ""
    
    if not hasattr(X, "col_names") or not hasattr(X, "row_names"):
        DIAGNOSIS = "No annotations."
        return False
    x = [x.upper() for x in X.row_names()]
    if sorted(x) != sorted(["ACCESSION", "DESCRIPTION", "CALL"]):
        DIAGNOSIS = "Improper row headers."
        return False
    x = [x.upper() for x in X.col_names()]
    if sorted(x) != sorted([tdf.SAMPLE_NAME, "DESCRIPTION", "SCALE_FACTOR"]):
        DIAGNOSIS = "Improper column headers."
        return False
    return True

def read(handle, hrows=None, hcols=None, datatype=float):
    from genomicode import filelib
    from genomicode import jmath
    from genomicode import Matrix
    import tab_delimited_format as tdf
    import const

    handle = filelib.openfh(handle)
    # Can't use iolib.split_tdf here because it does not handle empty
    # lines properly (which can occur if there is a file with no
    # samples).
    #data = iolib.split_tdf(handle.read())
    data = [x.rstrip("\r\n").split("\t") for x in handle]
    assert len(data) >= 3, "Invalid RES file."

    # Do some checking on the format.
    assert len(data[0]) == len(data[1])+1
    x = sorted([x.upper() for x in data[0][:2]])
    assert x == ["ACCESSION", "DESCRIPTION"]
    assert len(data[2]) == 1, "%d: %s" % (len(data[2]), repr(data[2]))

    # Parse out the number of genes and delete the row.
    num_genes = int(data[2][0])
    del data[2]
    assert len(data) == num_genes+2  # data + 2 headers

    # GenePattern creates files where the last column is all blank.
    # If this is the case, then delete it.
    blank_last_col = True
    x = [x[-1] for x in data if x[-1]]
    if not x:
        # Last column is all blank so delete it.
        data = [x[:-1] for x in data]

    # Parse the names of the samples.
    sample_names = []
    for i, x in enumerate(data[0][2:]):
        if i % 2:
            assert not x
        else:
            assert x
            sample_names.append(x)

    # Parse out the sample_description.
    sample_description = []
    for i, x in enumerate(data[1]):
        if i%2 == 0:
            assert not x
        else:
            assert x
            sample_description.append(x)
    assert len(sample_description) == len(sample_names)

    # Pull the scale factors out of the sample_description.
    # Some of the descriptions can be missing scale factors.
    scale_factors = [""] * len(sample_description)
    for i in range(len(sample_description)):
        x = sample_description[i]
        sf = "scale factor"
        j = x.lower().find(sf)
        if j < 0:
            continue
        assert x[j-1] == "/"
        assert x[j+len(sf)] == "="
        scale_factors[i] = float(sample_description[i][j+len(sf)+1:])
        sample_description[i] = sample_description[i][:j-1]

    # Parse out the description and accession columns.
    accession_header = data[0][0]
    description_header = data[0][1]
    accession = [x[0] for x in data[2:]]
    description = [x[1] for x in data[2:]]
    x = [x.upper() for x in data[0][:2]]
    if x == ["DESCRIPTION", "ACCESSION"]:
        accession_header, description_header = \
                          description_header, accession_header
        accession, description = description, accession
    assert (accession_header.upper(), description_header.upper()) == \
           ("ACCESSION", "DESCRIPTION")
    
    # Accession should be unique.
    x = {}.fromkeys(accession).keys()
    assert len(x) == len(accession)

    # Parse out the matrix and calls.
    matrix = []
    calls = []
    for row in data[2:]:
        row = row[2:]
        x0 = [x for (i, x) in enumerate(row) if i%2==0]
        x1 = [x for (i, x) in enumerate(row) if i%2==1]
        assert len(x0) == len(x1)
        for x in x1:
            assert x.upper() in ["A", "P", "M"], x
        matrix.append(x0)
        calls.append(x1)
    assert len(matrix) == num_genes
        
    # Should have some way of specifying no conversion.
    if datatype is None:
        convert_fn = None   # default
    elif datatype is int:
        convert_fn = jmath.safe_int
    elif datatype is float:
        convert_fn = jmath.safe_float
    else:
        convert_fn = datatype
        
    if convert_fn:
        matrix = [map(convert_fn, x) for x in matrix]

    row_names = {}
    col_names = {}
    row_order = data[0][:2] + ["CALL"]
    col_order = [tdf.SAMPLE_NAME, "DESCRIPTION", "SCALE_FACTOR"]
    
    row_names[accession_header] = accession
    row_names[description_header] = description
    # Store the calls as row annotations.  The gene annotation "CALL"
    # is a string of A, P, or M, with one call per sample.
    row_names["CALL"] = ["".join(x) for x in calls]
    
    col_names[tdf.SAMPLE_NAME] = sample_names
    col_names["DESCRIPTION"] = sample_description
    col_names["SCALE_FACTOR"] = scale_factors

    synonyms = {}
    synonyms[const.COL_ID] = tdf.SAMPLE_NAME
    synonyms[const.ROW_ID] = accession_header

    X = Matrix.InMemoryMatrix(
        matrix, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #X = Matrix.add_synonyms(X, synonyms)
    #is_matrix(X); print DIAGNOSIS
    assert is_matrix(X)
    return X

def write(X, handle):
    raise NotImplementedError
