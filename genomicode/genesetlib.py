"""

read_gmx
read_gmt
read_genesets   Read the genesets in a geneset file.
read_genes      Read a list of genes that belong in a geneset.

detect_format

"""

def read_gmx(filename, preserve_spaces=False):
    # yield name, description, list of genes
    import StringIO
    import filelib

    # Transpose this file and parse as gmt.
    matrix = [x for x in filelib.read_cols(filename)]
    assert len(matrix) >= 2
    t_matrix = []
    for j in range(len(matrix[0])):
        x = []
        for i in range(len(matrix)):
            # These lines may not be the same length as the names,
            # e.g. if the gene sets at the end have very few genes.
            if j >= len(matrix[i]):
                continue
            x.append(matrix[i][j])
        assert len(x) >= 2
        x1 = x[:2]
        x2 = [x.strip() for x in x[2:]]
        x = x1 + x2
        t_matrix.append(x)

    handle = StringIO.StringIO()
    for x in t_matrix:
        print >>handle, "\t".join(x)

    handle.seek(0)
    return read_gmt(handle, preserve_spaces=preserve_spaces)

## def read_gmx(filename):
##     # yield name, description, list of genes
##     import filelib
##     genesets = {}  # name -> list of genes
##     name2desc = {} # name -> description

##     iter = filelib.read_cols(filename)

##     # Read the names line.
##     names = iter.next()
##     for name in names:
##         genesets[name] = []

##     # Read the description line.
##     descs = iter.next()
##     assert len(descs) == len(names)
##     for name, desc in zip(names, descs):
##         name2desc[name] = desc

##     # Read each of the genes.
##     # These lines may not be the same length as the names, e.g. if the
##     # gene sets at the end have very few genes.
##     for cols in iter:
##         for i in range(len(cols)):
##             gene = cols[i].strip()
##             if not gene:
##                 continue
##             name = names[i]
##             genesets[name].append(gene)

##     for name in names:
##         desc = name2desc[name]
##         genes = genesets[name]
##         yield name, desc, genes

def read_gmt(filename, preserve_spaces=False):
    # yield name, description, list of genes
    import filelib
    import iolib

    # <name> <desc> [<gene>, <gene>, ...]
    x = filelib.openfh(filename).read()
    for cols in iolib.split_tdf(x):
        assert len(cols) >= 2
        name, description = cols[:2]
        x = cols[2:]
        x = [x.strip() for x in x]
        if not preserve_spaces:
            x = [x for x in x if x]
        genes = x
        #print preserve_spaces, len(genes), len(cols)
        yield name, description, genes

def read_genesets(filename):
    # yield name, description, list of genes
    if filename is None:
        return
    
    # Figure out the format of the geneset file.
    fmt = detect_format(filename)
    assert fmt is not None, "I could not figure out the format of: %s" % \
           filename
    if fmt == "GMX":
        read_fn = read_gmx
    elif fmt == "GMT":
        read_fn = read_gmt
    else:
        raise AssertionError, "Unknown format: %s" % fmt

    # Read the geneset file.
    for x in read_fn(filename):
        yield x

def read_genes(filename, *genesets):
    """Return a list of genes from the desired gene sets.  If no
    genesets are provided, then I will return a list of all the genes
    in the file.

    """
    if filename is None:
        return []
    
    # Read the geneset file.
    geneset2genes = {}
    all_genesets = []
    for x in read_genesets(filename):
        name, description, genes = x
        geneset2genes[name] = genes
        all_genesets.append(name)

    # Pull out the appropriate genes.
    if not genesets:
        genesets = all_genesets
    all_genes = []
    for gs in genesets:
        assert gs in geneset2genes, "Unknown geneset: %s" % gs
        genes = geneset2genes[gs]
        all_genes.extend(genes)

    # Get rid of duplicates while maintaining the original order of
    # the genes.
    i = 0
    seen = {}
    while i < len(all_genes):
        if all_genes[i] in seen:
            del all_genes[i]
        else:
            seen[all_genes[i]] = 1
            i += 1
    return all_genes

def _is_known_desc(desc):
    ldesc = desc.lower()
    if ldesc == "na":
        return True
    if ldesc.startswith("http"):
        return True
    if ldesc == "":
        return True
    return False

def detect_format(filename):
    # Return "GMX", "GMT", or None.
    import filelib
    import iolib

    x = filelib.openfh(filename).read()
    matrix = [x for x in iolib.split_tdf(x)]

    # GMX
    # <name1> <name2> ... <nameN>
    # <desc1> <desc2> ... <descN>
    # <gene>  <gene>  ... <gene>
    # <gene>   ""         <gene>
    #
    # GMT
    # <name1> <desc1> <gene> <gene>
    # <name2> <desc2> <gene>  ""
    # ...      ...     ...
    # <nameN> <descN> <gene> <gene>
    #
    # Examples of <desc>: na, http://

    # Clean up each of the cells of the matrix.
    for i in range(len(matrix)):
        matrix[i] = [x.strip() for x in matrix[i]]

    # Both the GMX or GMT formats may have different number of columns
    # in each row.  Align the rows so they have the same number of
    # columns.
    num_rows = [len(x) for x in matrix]
    if num_rows:
        assert min(num_rows) > 0
        max_rows = max(num_rows)
        for i in range(len(matrix)):
            x = [""]*(max_rows-len(matrix[i]))
            matrix[i] = matrix[i] + x

    # Cases:
    # 1.  Empty file               Broken
    # 2.  0 rows or 0 columns      Broken
    # 3.  1 row and 1 column       Broken
    # 4.  >= 2 rows, 1 column      GMX
    # 5.  1 row, >= 2 columns      GMT
    # 6.  2 rows and 2 columns     Check for known descriptions.
    # 7.  > 2 rows, > 2 columns    Check for spaces, then descriptions

    GMX, GMT = "GMX", "GMT"

    # Case 1.
    if not matrix:
        return None
    # Case 2.
    if not matrix[0]:
        return None
    nrow, ncol = len(matrix), len(matrix[0])
    # Case 3.
    if nrow == 1 and ncol == 1:
        return None
    # Case 4.
    if nrow >= 2 and ncol == 1:
        return GMX
    # Case 5.
    if nrow == 1 and ncol >= 2:
        return GMT
    # Case 6.
    if nrow == 2 and ncol == 2:
        desc01 = _is_known_desc(matrix[0][1])
        desc10 = _is_known_desc(matrix[1][0])
        if desc10 and not desc01:
            return GMX
        if desc01 and not desc10:
            return GMT
        return None
    # Case 7.
    assert nrow > 2 and ncol > 2

    # Check the extension of the file.
    if type(filename) is type(""):
        lfilename = filename.lower()
        if lfilename.endswith(".gmt"):
            return GMT
        if lfilename.endswith(".gmx"):
            return GMX

    # Check if there are any spaces in the first row or column.  If
    # the first row (col 3 onwards) contains any spaces, then it's not
    # a GMX file.  Vice versa for the first column.
    row1 = matrix[0]
    x = [x for x in row1[2:] if x == ""]
    spaces_in_row1 = (len(x) > 0)
    col1 = [x[0] for x in matrix]
    x = [x for x in col1[2:] if x == ""]
    spaces_in_col1 = (len(x) > 0)
    if spaces_in_row1 and spaces_in_col1:
        return None
    if not spaces_in_row1 and spaces_in_col1:
        return GMX
    if not spaces_in_col1 and spaces_in_row1:
        return GMT

    # Most likely, there will be different numbers of genes per gene
    # set.  This will lead to there being empty strings.  Check for a
    # pattern of empty strings and filled strings to distinguish GMX
    # or GMT.
    genes_left_aligned = True
    for i in range(nrow):
        is_left_aligned = True
        found_space = False
        row = matrix[i]
        for x in row[2:]:
            if x == "":
                found_space = True
            elif found_space:
                is_left_aligned = False
        if not is_left_aligned:
            genes_left_aligned = False
    genes_top_aligned = True
    for i in range(ncol):
        is_top_aligned = True
        found_space = False
        col = [x[i] for x in matrix]
        for x in col[2:]:
            if x == "":
                found_space = True
            elif found_space:
                is_top_aligned = False
        if not is_top_aligned:
            genes_top_aligned = False
    if not genes_top_aligned and not genes_left_aligned:
        return None
    if genes_top_aligned and not genes_left_aligned:
        return GMT
    if genes_left_aligned and not genes_top_aligned:
        return GMX

    # Check the descriptions to see if they match up.
    row2 = matrix[1]
    desc_row = True
    for x in row2:
        if not _is_known_desc(x):
            desc_row = False
    col2 = [x[1] for x in matrix]
    desc_col = True
    for x in col2:
        if not _is_known_desc(x):
            desc_col = False
    if desc_row and not desc_col:
        return GMX
    if desc_col and not desc_row:
        return GMT
    
    return None
        
                
def test_detect_format():
    from StringIO import StringIO
    
    TESTS = [
        ([], None),         # Empty file
        ([["foo"]], None),    # Too small.
        ([["name1", "desc1", "gene", "gene", "gene"],], "GMT"),
        ([["name1"],
          ["desc1"],
          ["gene"],
          ["gene"],
          ], "GMX"),
        ([["name1", "name2"],
          ["desc1", "desc2"],
          ], None),    # Same rows and columns.
        ([["name1", "name2"],
          ["na", "na"],
          ], "GMX"),
        ([["name1", "na"],
          ["name2", "na"],
          ], "GMT"),
        ([["name1", "name2", "name3"],
          ["desc1", "desc2", "desc3"],
          ["gene", "gene", "gene"],
          ], None),    # Same rows and columns.
        ([["name1", "name2", "name3"],
          ["na", "na", "na"],
          ["gene", "gene", "gene"],
          ], "GMX"),    # Same rows and columns.
        ([["name1", "name2", "name3"],
          ["desc1", "desc2", "desc3"],
          ["gene", "", "gene"],
          ], None), 
        ([["name1", "name2", "name3"],
          ["desc1", "desc2", "desc3"],
          ["", "", "gene"],
          ], "GMX"), 
        ([["name1", "name2", "name3"],
          ["desc1", "desc2", "desc3"],
          ["", "", "gene"],
          ], "GMX"), 
        ([["name1", "desc1", "gene"],
          ["name2", "desc2", "gene"],
          ["name3",],
          ], None), 
        ]

    header = "I", "GOLDEN", "DETECTED"
    print "\t".join(header)
    for i, x in enumerate(TESTS):
        #if i != 8:
        #    continue
        test, standard = x
        test_str = "\n".join(["\t".join(x) for x in test])+"\n"
        handle = StringIO(test_str)
        outcome = detect_format(handle)
        x = i, standard, outcome
        print "\t".join(map(str, x))


if __name__ == '__main__':
    test_detect_format()
