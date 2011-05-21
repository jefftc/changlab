"""Implement Output Description Format (ODF) from GenePattern.

Format described at:
http://www.broad.mit.edu/cancer/software/genepattern/tutorial/gp_fileformats.html

ODF 1.0
HeaderLines=7
Model= Dataset
DataLines= 3
COLUMN_TYPES: String String float float float *
COLUMN_DESCRIPTIONS: Sample from DFCI Sample from UK Sample from Children's
COLUMN_NAMES: Name Description Sample 1 Biopsy_2 Biopsy_4
RowNamesColumn=0
RowDescriptionsColumn=1

* First 5 lines required, but only first 2 lines must be in this
  order.
* Lines starting with # should be ignored.
* Rows are required to have a name and description.  The
  RowNamesColumn and RowDescriptionsColumn are optional, and it is
  unspecified what to do if these are not specified.  Appears to be
  0-based.
* Columns have a name, but the description is optional.  It is
  unspecified what to do if COLUMN_NAMES is not given.
* It is unspecified whether row names needs to be unique.
* COLUMN_NAMES gives the names for the row and description columns,
  while COLUMN_DESCRIPTIONS skips these.



Functions:
read
write
is_format
is_matrix

"""
import os

def is_format(locator_str):
    from genomicode import filefns
    if not filefns.exists(locator_str):
        return False

    # Read the first 2 non-comment lines.
    handle = filefns.openfh(locator_str)
    lines = []
    for line in handle:
        if line.startswith("#"):
            continue
        lines.append(line)
        if len(lines) >= 2:
            break
    handle.close()   # need to close it properly, or gunzip might not die.

    if len(lines) < 2:
        return False
    if not lines[0].startswith("ODF 1.0"):
        return False
    if not lines[1].startswith("HeaderLines"):
        return False
    return True

def is_matrix(X):
    if not hasattr(X, "col_headers") or not hasattr(X, "row_headers"):
        return False
    if len(X.row_headers()) != 2:
        return False
    col_headers = sorted(X.col_headers())
    if not col_headers:
        return False
    elif len(col_headers) == 1:
        if col_headers != ["TYPE"]:
            return False
    elif len(col_headers) == 2:
        if sorted(col_headers) != sorted(["TYPE", "DESCRIPTION"]):
            return False
    else:
        return False
    return True

def read(handle, datatype=None):
    from genomicode import filefns
    from genomicode import Matrix
    import const

    # datatype is not used here.  It is explicitly specified in the
    # format.
    handle = filefns.openfh(handle)

    # Read the header.  The format description doesn't specify whether
    # the names are case sensitive, so accept case insensitive names
    # by converting everything to uppercase.
    x = handle.readline().strip()
    assert x == "ODF 1.0", "Missing ODF version."
    x = handle.readline().strip().split("=")
    assert len(x) == 2
    assert x[0].upper() == "HEADERLINES"
    header_lines = int(x[1])
    assert header_lines >= 3 and header_lines <= 7, \
           "Invalid number of header lines."
    lines = [handle.readline() for i in range(header_lines)]
    lines = [x for x in lines if x]   # remove blank lines if EOF
    assert len(lines) == header_lines, "Wrong number of lines in header."

    # Parse the header lines.
    header = {}  # name -> value
    num_cols = None   # just the data, not the annotation headers.
    for line in lines:
        delimiter = "="
        if line.startswith("COLUMN"):
            delimiter = ":"
        assert delimiter in line, "Header missing delimiter '%s': %s" % (
            delimiter, line)
        name, value = line.split(delimiter)
        name, value = name.strip(" \r\n"), value.strip(" \r\n")
        if name.startswith("COLUMN"):
            value = value.split("\t")
            num_data = len(value)
            if name in ["COLUMN_TYPES", "COLUMN_NAMES"]:
                # Contains metadata describing the annotations.
                num_data = len(value)-2
            if num_cols is None:
                num_cols = num_data
            assert num_data == num_cols
        name = name.upper()
        header[name] = value
    header["DATALINES"] = int(header["DATALINES"])

    assert "MODEL" in header, 'Missing "Model" header.'
    assert "DATALINES" in header, 'Missing "DataLines" header.'
    assert "COLUMN_TYPES" in header, 'Missing "COLUMN_TYPES" header.'
    assert num_cols is not None  # Should come from COLUMN_TYPES.
    assert header["DATALINES"] >= 0 and header["DATALINES"] < 1E6

    # Read the data block.
    lines = [handle.readline() for i in range(header["DATALINES"])]
    lines = [x for x in lines if not x.startswith("#")]  # no comments
    data = [x.rstrip("\r\n").split("\t") for x in lines]
    # There might be leftover information in the file.  The format
    # does not describe how to handle this case.
    
    # Parse the column names out of the header.
    col_types = header["COLUMN_TYPES"]                    # required
    col_names = header.get("COLUMN_NAMES")                # optional
    col_descriptions = header.get("COLUMN_DESCRIPTIONS")  # optional

    # Make sure each line has the right number of columns.
    for x in data:
        assert len(x) == num_cols+2  # add the 2 columns of annotations
    # Convert the types of the data.
    for j in range(num_cols):
        coltype = col_types[j]
        ucoltype = coltype.upper()
        convert_fn = None
        if ucoltype == "STRING":
            pass
        elif ucoltype == "FLOAT":
            convert_fn = float
        else:
            raise AssertionError, "Unknown column type: %s" % coltype
        if not convert_fn:
            continue
        for i in range(len(data)):
            data[i][j] = convert_fn(data[i][j])

    # The first two columns are for the row names and description.
    col0 = [x[0] for x in data]
    col1 = [x[1] for x in data]
    head0, head1 = "Name", "Description"
    if col_names:
        head0, head1 = col_names[0], col_names[1]
    row_names, row_descriptions = col0, col1
    row_name_header, row_description_header = head0, head1
    if "ROWNAMESCOLUMN" in header:
        col = int(header["ROWNAMESCOLUMN"])
        assert col in [0, 1]
        if col == 1:
            row_names = col1
            row_name_header = head1
    if "ROWDESCRIPTIONSCOLUMN" in header:
        col = int(header["ROWDESCRIPTIONSCOLUMN"])
        assert col in [0, 1]
        if col == 0:
            row_descriptions = col0
            row_description_header = head0

    # Cut off the headers for the annotations.
    col_types = col_types[2:]
    if col_names:
        col_names = col_names[2:]

    matrix = [x[2:] for x in data]
    row_headers = [head0, head1]
    col_headers = None
    row_annots = {}
    col_annots = {}
    synonyms = {}

    row_annots[row_name_header] = row_names
    row_annots[row_description_header] = row_descriptions
    col_annots["TYPE"] = col_types
    if col_descriptions:
        col_annots["DESCRIPTION"] = col_descriptions
    synonyms[const.ROW_ID] = row_name_header
        
    X = Matrix.InMemoryMatrix(
        matrix, row_names=None, col_names=col_names,
        row_headers=row_headers, col_headers=col_headers,
        row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
    assert is_matrix(X)
    return X

def write(X, handle):
    import const
    
    print >>handle, "ODF 1.0"

    has_descriptions = "DESCRIPTIONS" in X.col_headers()
    has_names = not X.col_names()[0].startswith("C0")
    num_headers = 5 + int(has_descriptions) + int(has_names)

    print >>handle, "HeaderLines=%d" % num_headers
    print >>handle, "Model= Dataset"
    print >>handle, "DataLines= %s" % X.nrow()

    assert len(X.row_headers()) == 2
    x = ["String"]*2 + X.col_annots("TYPE")
    x = "\t".join(x)
    print >>handle, "COLUMN_TYPES: %s" % x

    if has_descriptions:
        x = "\t".join(X.col_annots("DESCRIPTION"))
        print >>handle, "COLUMN_DESCRIPTIONS: %s" % x
        
    if has_names:
        x = X.row_headers() + X.col_names()
        x = "\t".join(x)
        print >>handle, "COLUMN_NAMES: %s" % x

    col_names, col_desc = 0, 1
    name = X._synonyms[const.ROW_ID]
    if X.row_headers().index(name) != col_names:
        col_names, col_desc = col_desc, col_names
    assert X.row_headers().index(name) == col_names
    print >>handle, "RowNamesColumn=%d" % col_names
    print >>handle, "RowDescriptionsColumn=%d" % col_desc

    annots0 = X.row_annots(X.row_headers()[0])
    annots1 = X.row_annots(X.row_headers()[1])
    for i in range(X.nrow()):
        x = [annots0[i], annots1[i]] + X.value(i, None)
        print >>handle, "\t".join(map(str, x))
    handle.flush()
