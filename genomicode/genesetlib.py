"""

read_gmx        Each column is a gene set.
read_gmt        Each row is a gene set.
read_tdf
read_genesets   Read the genesets in a geneset file.
read_genes      Return a list of genes that belong in a specific geneset.

write_gmt
write_gmx

detect_format

score_geneset
score_geneset_I

clean_genes

"""
# Optimization: call _is_known_desc less.


DEBUG = None

class GeneSet:
    def __init__(self, name, description, genes):
        assert type(name) is type("")
        assert type(description) is type("") or description is None
        assert type(genes) is type([])

        if description is None or description.strip() == "":
            description = "na"
            
        # No duplicate genes, preserve order.
        nodup = []
        for g in genes:
            if g not in nodup:
                nodup.append(g)
        genes = nodup
        
        self.name = name
        self.description = description
        self.genes = genes


def write_gmt(filename, genesets):
    handle = filename
    if type(handle) is type(""):
        handle = open(filename, 'w')
    for gs in genesets:
        x = [gs.name, gs.description] + gs.genes
        print >>handle, "\t".join(map(str, x))


def write_gmx(filename, genesets):
    handle = filename
    if type(handle) is type(""):
        handle = open(filename, 'w')

    matrix = []
    for gs in genesets:
        x = [gs.name, gs.description] + gs.genes
        matrix.append(x)

    # Make sure each row of the matrix is the same length.
    maxlen = max([len(x) for x in matrix])
    for i in range(len(matrix)):
        x = [""] * (maxlen-len(matrix[i]))
        matrix[i].extend(x)

    for i in range(len(matrix[0])):
        x = [x[i] for x in matrix]
        print >>handle, "\t".join(map(str, x))


def read_gmx(filename, preserve_spaces=False, allow_duplicates=False):
    # yield name, description, list of genes
    import StringIO
    import filelib

    #matrix = [x for x in filelib.read_cols(filename)]
    matrix = filelib.read_all_cols(filename)
    assert len(matrix) >= 2
    
    # Transpose this file and parse as gmt.
    t_matrix = _transpose_gmx(matrix)
    handle = StringIO.StringIO()
    for x in t_matrix:
        print >>handle, "\t".join(x)
    handle.seek(0)
    
    return read_gmt(
        handle, preserve_spaces=preserve_spaces,
        allow_duplicates=allow_duplicates)


def _transpose_gmx(matrix):
    # Transpose a GMX (column-oriented) matrix to GMT (row-oriented).
    # The rows are not guaranteed to be the same length.
    
    # GMX format:
    # <gene set name>  ...
    # <description>    ...
    # <gene>           ...
    if not matrix:
        return []
    assert matrix[0]

    # Check if the rows are the same length.  If so, can use a faster
    # algorithm.
    all_same = True
    for i in range(1, len(matrix)):
        if len(matrix[0]) != len(matrix[i]):
            all_same = False
            break

    t_matrix = []
    if all_same:
        for j in range(len(matrix[0])):
            x = [row[j] for row in matrix]
            t_matrix.append(x)
    else:
        # Rows are not the same length.
        # Iterate over each of the columns (gene sets).
        for j in range(len(matrix[0])):
            x = []
            # Pull this column out of the matrix into variable x.
            for i in range(len(matrix)):
                # These lines may not be the same length as the names,
                # e.g. if the gene sets at the end have very few genes.
                if j >= len(matrix[i]):
                    continue
                x.append(matrix[i][j])
            assert len(x) >= 2, "Short column %d: %s" % (j, x)
            #x1 = x[:2]
            #x2 = [x.strip() for x in x[2:]]
            #x = x1 + x2
            t_matrix.append(x)

    #for i, x in enumerate(t_matrix):
    #    x = [x.strip() for x in x]
    #    t_matrix[i] = x
        
    return t_matrix


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


def read_gmt(filename, preserve_spaces=False, allow_duplicates=False):
    # yield name, description, list of genes
    # genes can be duplicated.
    import filelib
    #import iolib

    # <name> <desc> [<gene>, <gene>, ...]
    matrix = filelib.read_all_cols(filename)
    #matrix = [x for x in filelib.read_cols(filename)]
    #x = filelib.openfh(filename).read()
    #matrix = [x for x in iolib.split_tdf(x)]
    maxlen = max([len(x) for x in matrix]) - 2 # subtract name, description
    for cols in matrix:
        assert len(cols) >= 2
        name, description = cols[:2]
        x = cols[2:]
        x = [x.strip() for x in x]
        if not preserve_spaces:
            x = [x for x in x if x]
        genes = x
        
        if not allow_duplicates:
            # Remove the duplicates while preserving the order.
            nodup = []
            seen = {}
            for x in genes:
                if x in seen:
                    continue
                seen[x] = 1
                nodup.append(x)
            genes = nodup
        
        if preserve_spaces and len(genes) < maxlen:
            genes = genes + [""]*(maxlen-len(genes))

        #print preserve_spaces, len(genes), len(cols)
        yield name, description, genes


def read_tdf(filename, preserve_spaces=False, allow_duplicates=False,
             delimiter="\t", ignore_lines_startswith=None,
             yield_lines_startswith=None, nrows=None, strip=True):
    # yield name, description (always ""), list of genes.
    # preserve_spaces determines whether to remove blank annotations.
    # allow_duplicates determines whether to remove duplicate
    # annotations.  If ignore_lines_startswith is not None, then will
    # drop all lines that start with this substring.  If
    # yield_lines_startswith is not None, then will yield these lines
    # as a single string (not a tuple).
    import filelib

    matrix = filelib.read_all_cols(filename, delimiter=delimiter, nrows=nrows)
    #matrix = [x for x in
    #          filelib.read_cols(filename, delimiter=delimiter, nrows=nrows)]
    if strip:
        for i, x in enumerate(matrix):
            x = [x.strip() for x in x]
            matrix[i] = x
    if ignore_lines_startswith:
        # BUG: Will not work if this contains the delimiter.
        assert ignore_lines_startswith.find(delimiter) < 0
        matrix = [
            x for x in matrix if not x[0].startswith(ignore_lines_startswith)]
    if yield_lines_startswith:
        assert yield_lines_startswith.find(delimiter) < 0
        to_yield = [
            x for x in matrix if x[0].startswith(yield_lines_startswith)]
        to_yield = [delimiter.join(x) for x in to_yield]
        matrix = [
            x for x in matrix if not x[0].startswith(yield_lines_startswith)]
        for x in to_yield:
            yield x
    
    # Make sure the header row is as long as the longest row.
    # Otherwise, _transpose_gmx will cut off the columns at the end.
    x = [len(x) for x in matrix]
    if len(x) >= 2:
        header_len = x[0]
        row_len = max(x[1:])
        if header_len < row_len:
            x = [""] * (row_len-header_len)
            matrix[0] = matrix[0] + x

    # If matrix is empty, then nothing to yield.
    if not matrix:
        return
    # If there is only one row, then yield the columns.
    if len(matrix) == 1:
        for name in matrix[0]:
            yield name, "", []
        return

    # This will fail if there is not at least 2 rows.
    t_matrix = _transpose_gmx(matrix)
    if not matrix:
        return

    # Make sure the matrix is not empty.
    maxlen = max([len(x) for x in t_matrix]) - 1
    assert maxlen >= 0

    for i in range(len(t_matrix)):
        x = t_matrix[i]
        name, description, genes = x[0], "", x[1:]

        if not allow_duplicates:
            # Remove the duplicates while preserving the order.
            nodup = []
            seen = {}
            for x in genes:
                if x in seen:
                    continue
                seen[x] = 1
                nodup.append(x)
            genes = nodup
        
        if len(genes) < maxlen:
            genes = genes + [""]*(maxlen-len(genes))
        yield name, description, genes

    # BUG: This will fail is there is a geneset with no genes in it,
    # because read_gmx requires at least the name and description
    # lines.
    # BUG: Does not preserve spaces.
    #for x in read_gmx(filename, preserve_spaces=preserve_spaces):
    #    name, description, genes = x
    #    genes = [description] + genes
    #    description = ""
    #    yield name, description, genes


def read_genesets(
    filename, allow_tdf=False, allow_duplicates=False, preserve_spaces=False):
    # yield name, description, list of genes
    # If allow_tdf is True, will also parse filename if it is a
    # tab-delimited file where each column is a geneset.
    if filename is None:
        return
    
    # Figure out the format of the geneset file.
    fmt = detect_format(filename)
    if fmt == "GMX":
        read_fn = read_gmx
    elif fmt == "GMT":
        read_fn = read_gmt
    elif fmt:
        raise AssertionError, "Unknown format: %s" % fmt
    elif allow_tdf:
        read_fn = read_tdf
    else:
        raise AssertionError, \
              "I could not figure out the format of geneset file: %s" % \
           filename

    # Read the geneset file.
    for x in read_fn(
        filename, allow_duplicates=allow_duplicates,
        preserve_spaces=preserve_spaces):
        yield x


def read_genes(filename, *genesets, **keywds):
    """Return a list of genes from the desired gene sets.  If no
    genesets are provided, then I will return a list of all the genes
    in the file.

    """
    allow_tdf = keywds.get("allow_tdf", False)
    allow_duplicates = keywds.get("allow_duplicates", False)
    for k in keywds:
        assert k in ["allow_tdf", "allow_duplicates"]
    
    if filename is None:
        return []
    
    # Read the geneset file.
    geneset2genes = {}
    all_genesets = []
    for x in read_genesets(
        filename, allow_tdf=allow_tdf, allow_duplicates=allow_duplicates):
        name, description, genes = x
        geneset2genes[name] = genes
        all_genesets.append(name)

    # Pull out the appropriate genes.
    if not genesets:
        genesets = all_genesets
    all_genes = []
    for gs in genesets:
        if gs is None:
            continue
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
    global DEBUG
    import filelib
    #import iolib

    #x = filelib.openfh(filename).read()
    #matrix = [x for x in iolib.split_tdf(x)]
    matrix = filelib.read_all_cols(filename)
    #matrix = [x for x in filelib.read_cols(filename)]

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

    # If the last column is blank, remove it.
    last_column_blank = True
    for x in matrix:
        if x[-1] != "":
            last_column_blank = False
            break
    if last_column_blank:
        for i in range(len(matrix)):
            matrix[i] = matrix[i][:-1]

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
    # 7.  >= 2 rows, >= 2 columns  Check for spaces, then descriptions

    GMX, GMT = "GMX", "GMT"

    # Case 1.
    if not matrix:
        DEBUG = "detect_format case 1"
        return None
    # Case 2.
    if not matrix[0]:
        DEBUG = "detect_format case 2"
        return None
    nrow, ncol = len(matrix), len(matrix[0])
    # Case 3.
    if nrow == 1 and ncol == 1:
        DEBUG = "detect_format case 3"
        return None
    # Case 4.
    if nrow >= 2 and ncol == 1:
        DEBUG = "detect_format case 4"
        if _is_known_desc(matrix[1][0]):
            return GMX
        return None
    # Case 5.
    if nrow == 1 and ncol >= 2:
        DEBUG = "detect_format case 5"
        if _is_known_desc(matrix[0][1]):
            return GMT
        return None
    # Case 6.
    if nrow == 2 and ncol == 2:
        DEBUG = "detect_format case 6"
        desc01 = _is_known_desc(matrix[0][1])
        desc10 = _is_known_desc(matrix[1][0])
        if desc10 and not desc01:
            return GMX
        if desc01 and not desc10:
            return GMT
        return None
    # Case 7.
    DEBUG = "detect_format case 7"
    assert nrow >= 2 and ncol >= 2, "%d %d" % (nrow, ncol)

    # Check the extension of the file.
    if type(filename) is type(""):
        lfilename = filename.lower()
        if lfilename.endswith(".gmt"):
            DEBUG = "detect_format filename"
            return GMT
        if lfilename.endswith(".gmx"):
            DEBUG = "detect_format filename"
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
        DEBUG = "detect_format spaces in row1 and col1"
        return None
    if not spaces_in_row1 and spaces_in_col1:
        DEBUG = "detect_format spaces col"
        return GMX
    if not spaces_in_col1 and spaces_in_row1:
        DEBUG = "detect_format spaces row"
        return GMT

    # Most likely, there will be different numbers of genes per gene
    # set.  This will lead to there being empty strings.  Check for a
    # pattern of empty strings and filled strings to distinguish GMX
    # or GMT.

    # If there are spaces in the middle of the row, then the genes are
    # not left aligned.  E.g. gene <space> gene.
    genes_left_aligned = True
    for i in range(nrow):
        row = matrix[i][2:]
        I_space = [i for (i, x) in enumerate(row) if not x.strip()]
        I_char = [i for (i, x) in enumerate(row) if x.strip()]
        if not I_space or not I_char:
            continue
        if I_space[0] < I_char[-1]:
            genes_left_aligned = False
            break
    # Same, but with columns.
    genes_top_aligned = True
    for i in range(ncol):
        col = [x[i] for x in matrix[2:]]
        I_space = [i for (i, x) in enumerate(col) if not x.strip()]
        I_char = [i for (i, x) in enumerate(col) if x.strip()]
        if not I_space or not I_char:
            continue
        if I_space[0] < I_char[-1]:
            genes_top_aligned = False
            break
    if not genes_top_aligned and not genes_left_aligned:
        DEBUG = "detect_format alignment"
        return None
    if genes_top_aligned and not genes_left_aligned:
        DEBUG = "detect_format alignment (top)"
        return GMX
    if genes_left_aligned and not genes_top_aligned:
        DEBUG = "detect_format alignment (top left)"
        return GMT

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
        DEBUG = "detect_format descriptions"
        return GMX
    if desc_col and not desc_row:
        DEBUG = "detect_format descriptions"
        return GMT
    
    DEBUG = "detect_format default"
    return None


def score_geneset(MATRIX, pos_genes, neg_genes):
    # Return MATRIX_p, MATRIX_n, num_matches, list of scores (one for
    # each column), list of p-values.
    
    import matrixlib

    all_genes = pos_genes + neg_genes
    x = matrixlib.find_best_row_header(MATRIX, all_genes)
    header, num_found, found, not_found = x
    assert num_found, "I could not find any genes in the gene set."

    (pos_I, x) = MATRIX._index(row=pos_genes, row_header=header)
    (neg_I, x) = MATRIX._index(row=neg_genes, row_header=header)
    MATRIX_p = MATRIX.matrix(pos_I, None)
    MATRIX_n = MATRIX.matrix(neg_I, None)

    x = score_geneset_I(MATRIX._X, pos_I, neg_I)
    x, x, num_rows, scores, scores_o, pvalues = x
    return MATRIX_p, MATRIX_n, num_rows, scores, scores_o, pvalues


def wilcox_test(x, y):
    import jmath
    from jmath import R_fn, R_equals, R_var

    R = jmath.start_R()
    R_equals(x, "x")
    R_equals(y, "y")
    wt = R_fn("wilcox.test", R_var("x"), R_var("y"), RETVAL="wt")
    p_value = R("wt$p.value")[0]
    return p_value


def wilcox_test_each_col(X, Y):
    # X and Y are matrices with same numbers of columns (but might
    # have different rows).  Run the wilcoxon text across each of the
    # columns.  Return a list of p-values.
    import jmath
    from jmath import R_fn, R_equals, R_var

    assert X and Y
    assert X[0] and Y[0]
    assert len(X[0]) == len(Y[0])
    ncol = len(X[0])

    R = jmath.start_R()
    R_equals(X, "X")
    R_equals(Y, "Y")
    R("p.values <- rep(NA, %d)" % ncol)
    lines =  []
    lines.append("for(i in 1:%d) {" % ncol)
    lines.append("x <- wilcox.test(X[,i], Y[,i]);")
    lines.append("p.values[i] <- x$p.value;")
    lines.append("}")
    x = "\n".join(lines)
    R(x)
    #for i in range(ncol):
    #R("x <- wilcox.test(X[,%d], Y[,%d])" % (i+1, i+1))
    #R("p.values[%d] <- x$p.value" % (i+1))
    #R("}")
    p_values = list(R("p.values"))
    assert len(p_values) == ncol
    return p_values
    

def score_geneset_I(X, I_pos, I_neg):
    # Return X_p, X_n, num_matches, list of scores
    import jmath

    assert len(X)
    for i in range(len(X)):
        assert len(X[i]) == len(X[0])

    # Make a matrix with the positive genes.
    X_p = [X[i] for i in I_pos]

    # Make a matrix with the negative genes.
    X_n_orig = [X[i] for i in I_neg]
    X_n = []
    for x in X_n_orig:
        x = [-x for x in x]
        X_n.append(x)

    # Make a matrix with all other genes.
    I_other = [i for i in range(len(X)) if i not in I_pos and i not in I_neg]
    X_o = [X[i] for i in I_other]

    # Combine the positive and negative genes.
    X_pn = X_p + X_n

    # Calculate the scores.
    assert X_pn, "No genes in gene set"
    assert X_o, "No genes not in gene set"
    num_rows = len(X_pn)
    scores_pn = jmath.mean(X_pn, byrow=False)
    scores_o = jmath.mean(X_o, byrow=False)
    #scores_delta = [scores_pn[i]-scores_o[i] for i in range(len(scores_pn))]

    # DEBUGGING
    #scores_X = jmath.mean(X, byrow=False)
    #for i in range(len(scores_pn)):
    #    print "SAMP %03d All: %g (%d genes)" % (i, scores_X[i], len(X))
    #    print "SAMP %03d Geneset: %g (%d genes)" % (
    #        i, scores_pn[i], len(X_pn))
    #    print "SAMP %03d Not geneset: %g (%d genes)" % (
    #        i, scores_o[i], len(X_o))

    # Calculate the p-values.
    # This step is very slow due to the R call.
    #nc = len(X[0])
    #pvalues = [None] * nc
    #for j in range(nc):
    #    col_pn = [x[j] for x in X_pn]
    #    col_o = [x[j] for x in X_o]
    #    # Compare col_pn and col_o with wilcoxon test.
    #    p = wilcox_test(col_o, col_pn)
    #    pvalues[j] = p
    pvalues = wilcox_test_each_col(X_o, X_pn)
    
    return X_p, X_n_orig, num_rows, scores_pn, scores_o, pvalues


def int_if_possible(x):
    try:
        x = int(x)
    except ValueError:
        pass
    return x
        

def clean_genes(genes, delim=None):
    x = genes
    if delim is not None:
        x2 = []
        for x in x:
            x2.extend(x.split(delim))
        x = x2
    x = [x.strip() for x in x if x.strip()]
    x = [x for x in x if x != "---"]
    # Sort numerically, of possible.
    x = [int_if_possible(x) for x in x]
    x = sorted(x)
    # Only unique genes.
    seen = {}
    i = 0
    while i < len(x):
        if x[i] in seen:
            del x[i]
        else:
            seen[x[i]] = 1
            i += 1
    x = [str(x) for x in x]
    return x


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
