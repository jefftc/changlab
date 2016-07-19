#!/usr/bin/env python

# Functions:
# read_matrices
# has_missing_values
# num_missing_values
# assert_no_missing_values
# transpose_matrix
# transpose_nonmatrix
# correlate_matrix
# correlate_against_matrix
# correlate_some_vs_all
# calc_mean
# calc_sd
# calc_range
# calc_pca
# average_row_indexes
# group_expression_by_samplename
#
# parse_indexes
# parse_names
# parse_geneset
# parse_correlations
# _parse_file_gs
# _parse_file_string
# _parse_file_num
# _parse_file_annot
# _parse_file_num_annot
# _read_annot_file
#
# select_col_indexes
# select_col_ids
# select_col_genesets
# select_col_annotation
# select_col_regex
# select_col_random
# rename_col_id
# replace_col_ids
# relabel_col_ids
# append_col_ids
# add_col_ids
#
# reorder_col_indexes
# reorder_col_cluster
# reorder_col_alphabetical
# reorder_col_byfile
# remove_duplicate_cols
# remove_unnamed_cols
# rename_duplicate_cols
# toupper_col_ids
# hash_col_ids
# apply_re_col_ids
# add_prefix_col_ids
# add_suffix_col_ids
#
# tcga_normal_only
# tcga_primary_tumor_only
# tcga_metastasis_only
# tcga_relabel_patient_barcodes
# tcga_label_by_tissue_type
#
# select_row_indexes
# select_row_ids
# select_row_string
# select_row_numeric
# select_row_random
# select_row_genesets
# select_row_annotation
# select_row_numeric_annotation
# select_row_nonempty
# select_row_maxvalue
# filter_row_minvalue
# filter_row_maxvalue
# select_row_mean_value
# select_row_mean_var
# select_row_missing_values
# select_row_var
# select_row_delta
# select_row_fc
# select_row_num_samples_fc
# dedup_row_by_var
# reverse_rows
# reorder_row_indexes
# reorder_row_cluster
# reorder_row_cor
# reorder_row_mean
# rename_duplicate_rows
#
# align_rows
# align_cols
#
# add_row_id
# add_row_annot
# remove_row_annot
# rename_row_annot
# move_row_annot
# concat_row_annot
#
# rlog_blind
# calc_cpm
# set_min_value
# normalize_rows_to
# normalize_rows_sub
# center_genes_mean
# center_genes_median
# normalize_genes_var
# count_missing_values
# median_fill_genes
# zero_fill_genes
# impute_missing_values_knn
# add_missing_values
# loess_normalize
#
# _match_rownames_to_geneset       DEPRECATED
# _match_colnames_to_geneset       DEPRECATED
# _best_match_colnames_to_geneset  DEPRECATED
# _num_matches_to_geneset          DEPRECATED
# _align_geneset_to_matrix         DEPRECATED
# _align_matrix_to_geneset         DEPRECATED
# _find_col_header                 DEPRECATED
#
# _clean
# _intersect_indexes
# _dedup_indexes
# _run_forked



def _clean(s):
    s = s.replace("\t", " ")
    s = s.strip()
    if s.startswith('"') and s.endswith('"'):
        s = s[1:-1]
    return s


def read_matrices(filenames, skip_lines, read_as_csv, remove_comments,
                  clean_only, num_header_cols):
    import os
    import sys
    import csv
    import tempfile
    import arrayio
    from genomicode import filelib

    temp_files = []
    try:
        if read_as_csv or remove_comments or skip_lines:
            delimiter = "\t"
            if read_as_csv:
                delimiter = ","
            if remove_comments:
                assert delimiter not in remove_comments

            # Make a temporary files for each matrix file.
            for filename in filenames:
                x, file_ = tempfile.mkstemp(dir=".")
                os.close(x)
                temp_files.append(file_)
            assert len(filenames) == len(temp_files)

            # Clean up the files.
            for (infile, outfile) in zip(filenames, temp_files):
                inhandle = filelib.openfh(infile)
                if skip_lines:
                    skip_lines = int(skip_lines)
                    for x in range(skip_lines):
                        inhandle.readline()

                if read_as_csv:
                    inhandle = csv.reader(inhandle)
                else:
                    inhandle = filelib.read_cols(inhandle, delimiter=delimiter)

                outhandle = open(outfile, 'w')
                for cols in inhandle:
                    if not cols:
                        continue
                    if remove_comments and cols[0].startswith(remove_comments):
                        continue
                    # Clean up the data.
                    cols = [_clean(x) for x in cols]

                    outhandle.write("\t".join(cols) + "\n")
                outhandle.close()

            # Use the cleaned files.
            filenames = temp_files

        if clean_only:
            if len(filenames) != 1:
                raise NotImplementedError("Can not clean and merge.")
            sys.stdout.write(open(filenames[0]).read())
            sys.exit(0)

        # Read the files.
        matrices = []
        for filename in filenames:
            assert os.path.exists(filename), \
                "I could not find the file: %s" % filename
            fmt_module = arrayio.choose_format(filename)
            assert fmt_module, \
                "I could not figure out the format of file: %s" % filename
            x = fmt_module.read(filename, hcols=num_header_cols)
            matrices.append(x)
    finally:
        for f in temp_files:
            if os.path.exists(f):
                os.unlink(f)

    return fmt_module, matrices


def num_missing_values(MATRIX):
    num_missing = 0
    for x in MATRIX._X:
        for y in x:
            if y is None:
                num_missing += 1
    return num_missing


def has_missing_values(MATRIX):
    for x in MATRIX._X:
        if None in x:
            return True
    return False


def assert_no_missing_values(MATRIX):
    num = num_missing_values(MATRIX)
    assert num == 0, "Matrix has %d missing values." % num


def transpose_matrix(MATRIX, transpose):
    from genomicode import jmath
    from genomicode import Matrix
    from arrayio import const
    from arrayio import tab_delimited_format as tdf

    if transpose is None:
        return MATRIX

    x = transpose.split(",")
    assert len(x) == 2, "Format: <old row ID>,<new row ID>."
    old_row_id, new_row_id = x

    assert old_row_id in MATRIX.row_names(), "missing header : %s" % (
        old_row_id)
    #print MATRIX.col_names()
    #print MATRIX.row_names()
    #import sys; sys.exit(0)

    rnames = [x for x in MATRIX.row_names() if x != old_row_id]

    X = jmath.transpose(MATRIX._X)
    row_order = [new_row_id]
    col_order = [tdf.SAMPLE_NAME] + rnames
    row_names = {
        row_order[0] : MATRIX.col_names(const.COL_ID),
        }
    col_names = {
        col_order[0] : MATRIX.row_names(old_row_id),
        }
    for rn in rnames:
        col_names[rn] = MATRIX.row_names(rn)
    synonyms = {
        const.ROW_ID : row_order[0],
        const.COL_ID : col_order[0],
        }
    MATRIX_t = Matrix.InMemoryMatrix(
        X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_t


def transpose_nonmatrix(filename):
    # Transpose the contents of this tab-delimited file and print the
    # results.
    import os
    from genomicode import filelib
    from genomicode import iolib
    from genomicode import jmath
    from arrayio import tab_delimited_format as tdf

    assert os.path.exists(filename), "File not found: %s" % filename
    x = filelib.openfh(filename).read()
    x = iolib.split_tdf(x, strip=True)
    data = tdf._clean_tdf(x)
    # Make sure each line has the same number of columns.
    max_cols = max([len(x) for x in data])
    for i in range(len(data)):
        x = data[i]
        if len(x) >= max_cols:
            continue
        n = max_cols - len(x)
        x = x + [""]*n
        assert len(x) == max_cols
        data[i] = x
    data_t = jmath.transpose(data)
    for x in data_t:
        print "\t".join(x)


def correlate_matrix(MATRIX, correlate):
    from genomicode import jmath
    from genomicode import Matrix
    from arrayio import const

    if not correlate:
        return MATRIX

    X_cor = jmath.cor(MATRIX._X, byrow=0)
    assert len(X_cor) == MATRIX.ncol()

    row_order = [MATRIX.col_names()[0]]
    col_order = [MATRIX.col_names()[0]]
    row_names = {
        row_order[0] : MATRIX.col_names(const.COL_ID),
        }
    col_names = {
        col_order[0] : MATRIX.col_names(const.COL_ID),
        }
    synonyms = {
        const.ROW_ID : row_order[0],
        const.COL_ID : col_order[0],
        }
    MATRIX_cor = Matrix.InMemoryMatrix(
        X_cor, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_cor


def correlate_against_matrix(MATRIX, correlate):
    from genomicode import jmath
    from genomicode import Matrix
    from genomicode import hashlib
    from arrayio import const
    import arrayio

    if not correlate:
        return MATRIX

    MATRIX2 = arrayio.read(correlate)

    assert MATRIX.nrow() == MATRIX2.nrow()
    assert MATRIX.ncol() == MATRIX2.ncol()

    # Correlate each of the rows.
    cors = [""] * MATRIX.nrow()
    for i in range(MATRIX.nrow()):
        x = MATRIX._X[i]
        y = MATRIX2._X[i]
        # Ignore None.
        I1 = [j for (j, a) in enumerate(x) if a is not None]
        I2 = [j for (j, a) in enumerate(y) if a is not None]
        I = [j for j in I1 if j in I2]
        if not I:
            continue
        x = [x[j] for j in I]
        y = [y[j] for j in I]
        cors[i] = jmath.cor_list(x, y, safe=1)
    X_cor = []
    for x in cors:
        X_cor.append([x])

    x = MATRIX.row_names() + MATRIX2.row_names()
    x = hashlib.uniquify_by_num(x)
    
    row_order = x
    col_order = [MATRIX.col_names()[0]]
    # XXX need to be unique.
    row_names = {}
    for i, rn in enumerate(MATRIX.row_names()):
        row_names[row_order[i]] = MATRIX.row_names(rn)
    for i, rn in enumerate(MATRIX2.row_names()):
        i += len(MATRIX.row_names())
        row_names[row_order[i]] = MATRIX2.row_names(rn)
    col_names = {
        col_order[0] : ["Correlation"]
        }

    synonyms = {
        const.ROW_ID : row_order[0],
        const.COL_ID : col_order[0],
        }
    MATRIX_cor = Matrix.InMemoryMatrix(
        X_cor, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_cor


def correlate_some_vs_all(MATRIX, ids):
    import itertools
    from genomicode import jmath
    from genomicode import Matrix
    from arrayio import const
    from arrayio import tab_delimited_format as tdf

    if not ids:
        return MATRIX

    ids = [x.strip() for x in ids]
    ids = ",".join(ids)
    I = parse_names(MATRIX, True, ids)
    assert I, "I could not find any rows: %s" % ids

    # Just use the first column as the name.
    assert MATRIX._row_order
    HEADER = MATRIX._row_order[0]
    row_names = MATRIX.row_names(HEADER)

    jmath.start_R()
    data_tab = []
    for (i, j) in itertools.product(I, range(MATRIX.nrow())):
        row_i = MATRIX._X[i][:]
        row_j = MATRIX._X[j][:]
        assert len(row_i) == len(row_j), "%d %d %d %d" % (
            i, j, len(row_i), len(row_j))

        # Remove missing values.
        k = 0
        while k < len(row_i):
            if row_i[k] is None or row_j[k] is None:
                del row_i[k]
                del row_j[k]
            else:
                k += 1
        assert len(row_i) == len(row_j)
        if not row_i:  # all missing values
            continue

        mean_i = jmath.mean(row_i)
        mean_j = jmath.mean(row_j)

        jmath.R_equals(row_i, "X")
        jmath.R_equals(row_j, "Y")
        jmath.R_fn("cor.test", row_i, row_j, method="pearson", RETVAL="x")
        estimate = jmath.R("x$estimate")[0]
        p_value = jmath.R("x$p.value")[0]
        x = row_names[i], row_names[j], len(row_i), mean_i, mean_j, \
            estimate, p_value, 0.0
        data_tab.append(x)

    # Calculate the false discovery rate.
    p_values = [x[6] for x in data_tab]
    fdr = jmath.cmh_fdr_bh(p_values)
    for i in range(len(data_tab)):
        x = data_tab[i]
        x = list(x[:-1]) + [fdr[i]]
        data_tab[i] = x

    x1 = "%s 1" % HEADER
    x2 = "%s 2" % HEADER
    header = [
        x1, x2, "Length", "Mean 1", "Mean 2", "Estimate", "p value", "FDR"]

    s1 = [x[0] for x in data_tab]
    s2 = [x[1] for x in data_tab]
    X_cor = [x[2:] for x in data_tab]

    row_order = header[:2]
    col_order = [tdf.SAMPLE_NAME]
    row_names = {
        row_order[0] : s1,
        row_order[1] : s2,
        }
    col_names = {
        col_order[0] : header[2:]
        }
    synonyms = {
        const.ROW_ID : row_order[0],
        const.COL_ID : col_order[0],
        }
    MATRIX_cor = Matrix.InMemoryMatrix(
        X_cor, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_cor


def calc_mean(MATRIX, calc_mean):
    from genomicode import jmath
    
    if not calc_mean:
        return MATRIX

    mean = [None] * len(MATRIX._X)
    # Calculate the mean of each row.  Handle missing values.
    for i in range(len(MATRIX._X)):
        x = MATRIX._X[i]
        x = [x for x in x if x]
        mean[i] = jmath.mean(x)
    assert len(mean) == len(MATRIX._X)

    header = "Mean"
    i = 1
    while header in MATRIX._row_names:
        header = "Mean %d" % i
        i += 1
    MATRIX_new = MATRIX.matrix()
    MATRIX_new._row_order.append(header)
    MATRIX_new._row_names[header] = mean
    return MATRIX_new


def calc_sd(MATRIX, calc_sd):
    from genomicode import jmath
    
    if not calc_sd:
        return MATRIX

    sd = [None] * len(MATRIX._X)
    # Calculate the standard deviation of each row.  Handle missing values.
    for i in range(len(MATRIX._X)):
        x = MATRIX._X[i]
        x = [x for x in x if x]
        sd[i] = jmath.stddev(x)
    assert len(sd) == len(MATRIX._X)

    header = "Standard Deviation"
    i = 1
    while header in MATRIX._row_names:
        header = "Standard Deviation %d" % i
        i += 1
    MATRIX_new = MATRIX.matrix()
    MATRIX_new._row_order.append(header)
    MATRIX_new._row_names[header] = sd
    return MATRIX_new


def calc_range(MATRIX, calc_range):
    if not calc_range:
        return MATRIX

    delta = [None] * len(MATRIX._X)
    # Calculate the standard deviation of each row.  Handle missing values.
    for i in range(len(MATRIX._X)):
        x = MATRIX._X[i]
        x = [x for x in x if x]
        delta[i] = max(x) - min(x)
    assert len(delta) == len(MATRIX._X)

    header = "Range"
    i = 1
    while header in MATRIX._row_names:
        header = "Range %d" % i
        i += 1
    MATRIX_new = MATRIX.matrix()
    MATRIX_new._row_order.append(header)
    MATRIX_new._row_names[header] = delta
    return MATRIX_new


def calc_pca_rows(MATRIX, calc_pca_rows):
    if not calc_pca_rows:
        return MATRIX
    from arrayio import const, tdf
    from genomicode import Matrix
    from genomicode import pcalib
    from genomicode import jmath

    # Calculate the principal components.
    K = min(MATRIX.nrow(), MATRIX.ncol())
    pc, perc_var = pcalib.svd_project_cols(MATRIX._X, K)
    # Make pc a K x ncol matrix.
    pc = jmath.transpose(pc)
    assert len(pc) == K
    assert len(pc[0]) == MATRIX.ncol()
    assert len(perc_var) == K

    pc_name = ["PC %d" % i for i in range(1, K+1)]
    
    row_order = ["PC", "Perc Variance"]
    col_order = [tdf.SAMPLE_NAME]
    row_names = {
        row_order[0] : pc_name,
        row_order[1] : map(str, perc_var),
        }
    col_names = {
        col_order[0] : MATRIX.col_names(const.COL_ID),
        }
    synonyms = {
        const.ROW_ID : row_order[0],
        const.COL_ID : col_order[0],
        }
    MATRIX_pca = Matrix.InMemoryMatrix(
        pc, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_pca


def average_row_indexes(MATRIX, indexes, count_headers):
    if not indexes:
        return MATRIX
    from genomicode import Matrix
    from genomicode import jmath
    
    I = parse_indexes(MATRIX, True, indexes, count_headers)
    X = MATRIX.slice(I, None)
    average = jmath.mean(X, byrow=0)

    # Add average to the Matrix.
    row_order = MATRIX._row_order[:]
    col_order = MATRIX._col_order[:]
    row_names = MATRIX._row_names.copy()
    col_names = MATRIX._col_names.copy()
    synonyms = MATRIX._synonyms.copy()

    X = MATRIX._X[:]
    X.insert(0, average)
    for n, values in row_names.iteritems():
        values.insert(0, "AVERAGE")
    M = Matrix.InMemoryMatrix(
        X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return M
    

def group_expression_by_samplename(MATRIX, group):
    from genomicode import Matrix
    from arrayio import const
    from arrayio import tab_delimited_format as tdf

    if not group:
        return MATRIX

    # Pull out the expression values and sample names.
    I = parse_names(MATRIX, True, group)
    assert I, "I could not find the row: %s" % group
    assert len(I) <= 1, "Multiple row matches: %s" % group
    values = MATRIX.value(I[0], None)

    assert tdf.SAMPLE_NAME in MATRIX.col_names(), "Missing sample names"
    sample_names = MATRIX.col_names(tdf.SAMPLE_NAME)

    # Find the groups and assign the expression values to each group.
    groups = sorted({}.fromkeys(sample_names))
    assert len(groups) <= len(sample_names), "No common sample names."
    group2values = {}
    for name, value in zip(sample_names, values):
        if name not in group2values:
            group2values[name] = []
        group2values[name].append(value)
    x = [len(x) for x in group2values.itervalues()]
    max_group_len = max(x)

    # Make the matrix.
    X = [[None]*len(groups) for i in range(max_group_len)]
    for i in range(len(groups)):
        values = group2values[groups[i]]
        for j in range(len(values)):
            X[j][i] = values[j]

    row_order = []
    col_order = [tdf.SAMPLE_NAME]
    row_names = {
        }
    col_names = {
        col_order[0] : groups,
        }
    synonyms = {
        const.COL_ID : col_order[0],
        }

    MATRIX_g = Matrix.InMemoryMatrix(
        X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return MATRIX_g


def parse_indexes(MATRIX, is_row, indexes_str, count_headers):
    # is_row is a boolean for whether these are row indexes.  Takes
    # 1-based indexes and returns a list of 0-based indexes.
    #
    # Example inputs:
    # 5
    # 1,5,10
    # 1-99,215-300
    from genomicode import parselib

    max_index = MATRIX.nrow()
    num_headers = len(MATRIX._col_names)
    if not is_row:
        max_index = MATRIX.ncol()
        num_headers = len(MATRIX._row_names)

    I = []
    for s, e in parselib.parse_ranges(indexes_str):
        if count_headers:
            s, e = s - num_headers, e - num_headers
        assert s >= 1
        s, e = s - 1, min(e, max_index)
        I.extend(range(s, e))

    # Remove duplicated indexes.  Need to preserve order.
    nodup = []
    for i in I:
        if i not in nodup:
            nodup.append(i)
    I = nodup

    return I


def parse_names(MATRIX, is_row, s):
    # Examples:
    # E2F1,E2F3,PCNA
    names = s.split(",")
    params = {"row": names}
    if not is_row:
        params = {"col": names}

    I_row, I_col = MATRIX._index(**params)
    I = I_row
    if not is_row:
        I = I_col
    return I


def parse_geneset(MATRIX, is_row, geneset):
    # Return a list of indexes that match the desired gene sets.
    from genomicode import genesetlib

    if not geneset:
        return []
    filename, genesets = _parse_file_gs(geneset)
    assert len(genesets) >= 1, "I could not parse a gene set from: %s" % \
           geneset

    keywds = {"allow_tdf": True}
    genes = genesetlib.read_genes(filename, *genesets, **keywds)
    # No blank genes.
    genes = [x for x in genes if x.strip()]
    params = {"row": genes}
    if not is_row:
        params = {"col": genes}
    I_row, I_col = MATRIX._index(**params)
    I = I_row
    if not is_row:
        I = I_col
    return I


def parse_correlations(correlations):
    # correlations is a comma-separated list of numbers.  Return a
    # list of numbers.
    x = correlations.split(",")
    x = [float(x) for x in x]
    return x


def _parse_file_gs(geneset):
    # Parse a geneset specified by the user.  geneset is in the format
    # of <filename>[,<geneset>,<geneset>,...].  Return a tuple of
    # <filename>, list of <geneset> (or empty list).
    # XXX what happens if this is an empty list?
    x = geneset.split(",")
    assert len(x) >= 1
    filename, genesets = x[0], x[1:]
    return filename, genesets


def _parse_file_string(annotation):
    # annotations is in the format:
    # <header>,<value>[,<value,...]
    # Return a tuple of <header>, list <value>.
    x = annotation.split(",")
    assert len(x) >= 2
    header, values = x[0], x[1:]
    return header, values


def _parse_file_num(annotation):
    # annotations is in the format:
    # <header>,<value>[,<value,...]
    # <value> can optionally start with modifiers "<", "<=", ">=" or ">".
    # Return a tuple of <header>, list of (<modifier>, <value>).
    # If no modifier is specified, then <modifier> is "=".  <value> is
    # a floating point number.
    x = annotation.split(",")
    assert len(x) >= 2
    header, values = x[0], x[1:]
    for i in range(len(values)):
        value = values[i]
        modifier = "="
        if value[:2] in ["<=", ">="]:
            modifier = value[:2]
            value = value[2:]
        elif value[0] in "<>=":
            modifier = value[0]
            value = value[1:]
        value = float(value)
        values[i] = modifier, value
    return header, values


def _parse_file_annot(annotation):
    # annotations is in the format:
    # <txt_file>,<header>,<value>[,<value>,...]
    # Return a tuple of <filename>, <header>, list of <value>.
    x = annotation.split(",")
    assert len(x) >= 3, "Format: <txt_file>,<header>,<value>[,<value>,...]"
    filename, header, values = x[0], x[1], x[2:]
    return filename, header, values


def _read_annot_file(filename):
    # Return (header2annots, all_headers, all_annots).
    import os
    from genomicode import genesetlib

    assert os.path.exists(filename), "I could not find annotation file: %s" % \
        filename
    header2annots = {}
    all_headers = []
    all_annots = []
    for x in genesetlib.read_tdf(
            filename, preserve_spaces=True, allow_duplicates=True):
        name, x, annots = x
        header2annots[name] = annots
        all_headers.append(name)
        all_annots.append(annots)
        #if name == "Distance":
        #    print "Found"
        #    x = [x for x in annots if str(x) == ""]
        #    print "NUM BLANK", len(x)

    assert all_headers
    for annots in all_annots[1:]:
        assert len(annots) == len(all_annots[0])
    return header2annots, all_headers, all_annots


def _parse_file_num_annot(annotation):
    # annotations is in the format:
    # <txt_file>,<header>,<value>[,<value,...]
    # <value> can optionally start with modifiers "<", "<=", ">=" or ">".
    # Return a tuple of <filename>, <header>, list of (<modifier>, <value>).
    # If no modifier is specified, then <modifier> is "=".  <value> is
    # a floating point number.
    x = annotation.split(",")
    assert len(x) >= 3
    filename, header, values = x[0], x[1], x[2:]
    for i in range(len(values)):
        value = values[i]
        modifier = "="
        if value[:2] in ["<=", ">="]:
            modifier = value[:2]
            value = value[2:]
        elif value[0] in "<>=":
            modifier = value[0]
            value = value[1:]
        value = float(value)
        values[i] = modifier, value
    return filename, header, values


def select_col_indexes(MATRIX, indexes_list, count_headers):
    if not indexes_list:
        return None
    indexes = []
    for x in indexes_list:
        x = parse_indexes(MATRIX, False, x, count_headers)
        indexes.extend(x)
    return indexes


def select_col_ids(MATRIX, ids):
    # IDs is a list of gene IDs or names.  Combine them into a single string.
    ids = [x.strip() for x in ids]
    ids = ",".join(ids)
    if not ids:
        return None
    return parse_names(MATRIX, False, ids)


def select_col_genesets(MATRIX, genesets):
    if not genesets:
        return None
    I = []
    for geneset in genesets:
        x = parse_geneset(MATRIX, False, geneset)
        I.extend(x)
    # No duplicates.  Keep current order.
    i = 0
    while i < len(I):
        if I[i] in I[:i]:
            del I[i]
        else:
            i += 1
    return I



def select_col_annotation(MATRIX, col_annotation):
    # Format: list of [<txt_file>,<header>,<value>[,<value,...]]
    from genomicode import matrixlib

    if not col_annotation:
        return None

    I_all = []
    for i, col_annot in enumerate(col_annotation):
        x = _parse_file_annot(col_annot)
        filename, header, values = x

        # Read the annotations.
        x = _read_annot_file(filename)
        header2annots, all_headers, all_annots = x

        # Align the annotations to the matrix file.
        x = matrixlib.align_cols_to_many_annots(
            MATRIX, all_annots, hash=True, get_indexes=True)
        I_matrix, I_annots, index = x
        assert len(I_matrix) == len(I_annots)

        # Select only the columns where the annotation contains one of
        # these values.
        assert header in header2annots, "Missing header: %s" % header
        annots = header2annots[header]
        I = [im for (im, ia) in zip(I_matrix, I_annots)
             if annots[ia] in values]
        if i == 0:
            I_all = I
        else:
            I_all = sorted(set(I_all).intersection(I))

    return I_all


def select_col_regex(MATRIX, col_regex):
    import re
    from arrayio import tab_delimited_format as tdf

    if not col_regex:
        return None

    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    names = MATRIX.col_names(tdf.SAMPLE_NAME)
    I = []
    for i in range(len(names)):
        m = re.search(col_regex, names[i])
        if m:
            I.append(i)
    return I


def select_col_random(MATRIX, col_random):
    import random
    if col_random is None:
        return None
    assert col_random >= 1 and col_random < MATRIX.ncol()

    I = range(MATRIX.ncol())
    random.shuffle(I)
    I = sorted(I[:col_random])
    return I


def select_col_numeric_value(MATRIX, col_numeric):
    # Format: <row_id>,<value>[,<value,...]
    if not col_numeric:
        return None

    x = _parse_file_num(col_numeric)
    row_id, values = x

    I = parse_names(MATRIX, True, row_id)
    assert I, "I could not find row_id: %s" % row_id
    assert len(I) == 1, "Multiple matches for %s" % row_id

    signal = MATRIX._X[I[0]]
    I = []
    for i in range(len(signal)):
        x = float(signal[i])
        match = True
        for (modifier, value) in values:
            if modifier == "=":
                match = abs(x-value) < 1E-10
            elif modifier == "<":
                match = (x < value)
            elif modifier == ">":
                match = (x > value)
            elif modifier == "<=":
                match = (x <= value)
            elif modifier == ">=":
                match = (x >= value)
            else:
                raise AssertionError("Unknown modifier: %s" % modifier)
        # Accepts matches for any of the values.
        if match:
            I.append(i)
    return I


def rename_col_id(MATRIX, rename_list, ignore_missing):
    # rename_list is list of strings in format of: <from>,<to>.
    import arrayio

    if not rename_list:
        return MATRIX

    rename_all = []  # list of (from_str, to_str)
    for rename_str in rename_list:
        x = rename_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        from_str, to_str = x
        rename_all.append((from_str, to_str))

    name = arrayio.COL_ID
    if name not in MATRIX._col_names:
        name = MATRIX._synonyms[name]
    assert name in MATRIX._col_names, "I can not find the sample names."
    names = MATRIX.col_names(name)

    if not ignore_missing:
        missing = []
        for (from_str, to_str) in rename_all:
            if from_str not in names:
                missing.append(from_str)
        x1 = ", ".join(missing)
        x2 = names
        if len(x2) > 5:
            x2 = x2[:5] + ["..."]
        x2 = ", ".join(x2)
        assert not missing, "Can't find col ids (%s) in %s" % (x1, x2)

    MATRIX_new = MATRIX.matrix()
    for i in range(len(names)):
        for (from_str, to_str) in rename_all:
            if names[i] == from_str:
                names[i] = to_str
    MATRIX_new._col_names[name] = names

    return MATRIX_new


def replace_col_ids(MATRIX, replace_list):
    # replace_list is list of strings in format of: <from>,<to>.
    import arrayio

    if not replace_list:
        return MATRIX

    replace_all = []  # list of (from_str, to_str)
    for replace_str in replace_list:
        x = replace_str.split(",")
        assert len(x) == 2, "format should be: <from>,<to>"
        from_str, to_str = x
        replace_all.append((from_str, to_str))

    MATRIX_new = MATRIX.matrix()
    name = arrayio.COL_ID
    if name not in MATRIX_new._col_names:
        name = MATRIX_new._synonyms[name]
    assert name in MATRIX_new._col_names, "I can not find the sample names."
    x = MATRIX_new.col_names(name)
    for (from_str, to_str) in replace_all:
        x = [x.replace(from_str, to_str) for x in x]
    MATRIX_new._col_names[name] = x

    return MATRIX_new


def relabel_col_ids(MATRIX, geneset, ignore_missing):
    if not geneset:
        return MATRIX

    import os
    import sys
    import arrayio
    from genomicode import genesetlib
    from genomicode import matrixlib
    
    filename, genesets = _parse_file_gs(geneset)
    assert len(genesets) == 1, "Multiple genesets specified: %s." % geneset

    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    ext = os.path.splitext(filename)[1].lower()

    # Bug: sometimes will mis-identify TDF files as GMX.  The first
    # row will be interpreted as a description instead of a gene (or
    # annotation).  If the extension of the file isn't gmx or gmt,
    # then assume it's some sort of tdf file.  Also may mis-identify
    # as GMT (via detect_format alignment (top)).
    #if not genesetlib._is_known_desc(description) and \
    #       ext not in [".gmx", ".gmt"]:
    if ext not in [".gmx", ".gmt"]:
        read_fn = genesetlib.read_tdf
    else:
        fmt = genesetlib.detect_format(filename)
        if fmt == "GMX":
            read_fn = genesetlib.read_gmx
        elif fmt == "GMT":
            read_fn = genesetlib.read_gmt
        elif fmt:
            raise AssertionError, "Unknown format: %s" % fmt
        raise AssertionError, \
              "I could not figure out the format of geneset file: %s" % \
           filename

    for x in read_fn(
        filename, allow_duplicates=True, preserve_spaces=True):
        geneset, description, genes = x

        geneset2genes[geneset] = genes
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Make sure all the genes have the same length.  Otherwise,
    # something might be broken.
    assert all_genes
    for x in all_genes:
        assert len(x) == len(all_genes[0]), "genesets not aligned"

    # Find an alignment between the sample names and the genesets.
    x = matrixlib.align_cols_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True)
    I_matrix, I_geneset, index = x
    if len(I_matrix) == 0:
        raise AssertionError("Matrixes doesn't match any gene sets.")
    elif len(I_matrix) != MATRIX.ncol() and not ignore_missing:
        # No matches.  Try to diagnose.
        gs = all_genesets[index]
        nm = len(I_matrix)
        missing = []
        col_names = MATRIX.col_names(arrayio.COL_ID)
        missing = [col_names[i] for i in range(len(col_names))
                   if i not in I_matrix]
        print >>sys.stderr, \
            "Matrix best matches column '%s' [%d:%d]." % (
            gs, nm, MATRIX.ncol())
        MAX_TO_SHOW = 5
        if len(missing) > MAX_TO_SHOW:
            print >>sys.stderr, \
                "Missing (showing %d of %d):" % (MAX_TO_SHOW, len(missing))
        else:
            print >>sys.stderr, "Missing from annotation:"
        for x_ in missing[:MAX_TO_SHOW]:
            print >>sys.stderr, x_
        raise AssertionError("I could not match the matrix to a geneset.")

    # Add the new column names to the MATRIX.
    MATRIX_new = MATRIX.matrix()
    name = arrayio.COL_ID
    if name not in MATRIX_new._col_names:
        name = MATRIX_new._synonyms[name]
    assert name in MATRIX_new._col_names, "I can not find the sample names."
    names = MATRIX_new.col_names(name)
    gs = genesets[0]
    assert gs in geneset2genes, "Unknown header: %s" % gs
    genes = geneset2genes[gs]
    assert len(I_matrix) == len(I_geneset)
    assert max(I_geneset) < len(genes)
    assert max(I_matrix) < len(names)
    for i in range(len(I_matrix)):
        names[I_matrix[i]] = genes[I_geneset[i]]
    MATRIX_new._col_names[name] = names

    return MATRIX_new


def append_col_ids(MATRIX, geneset, ignore_missing):
    import os
    import sys
    import arrayio
    from genomicode import genesetlib
    from genomicode import matrixlib

    if not geneset:
        return MATRIX
    filename, genesets = _parse_file_gs(geneset)
    assert len(genesets) == 1

    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    ext = os.path.splitext(filename)[1].lower()
    for x in genesetlib.read_genesets(
            filename, allow_tdf=True, allow_duplicates=True):
        geneset, description, genes = x

        # Bug: sometimes will mis-identify TDF files as GMX.  The
        # first row will be interpreted as a description instead of a
        # gene (or annotation).  If the extension of the file isn't
        # gmx or gmt, and the description looks like a gene name, then
        # assume it's some sort of tdf file that has been
        # misidentified.
        #if not genesetlib._is_known_desc(description) and \
        #       ext not in [".gmx", ".gmt"]:
        if description and ext not in [".gmx", ".gmt"]:
            genes = [description] + genes

        geneset2genes[geneset] = genes
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Find an alignment between the sample names and the genesets.
    x = matrixlib.align_cols_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True)
    I_matrix, I_geneset, index = x
    if len(I_matrix) == 0:
        raise AssertionError("Matrixes doesn't match any gene sets.")
    elif len(I_matrix) != MATRIX.ncol() and not ignore_missing:
        # No matches.  Try to diagnose.
        gs = all_genesets[index]
        nm = len(I_matrix)
        missing = []
        col_names = MATRIX.col_names(arrayio.COL_ID)
        missing = [col_names[i] for i in range(len(col_names))
                   if i not in I_matrix]
        print >>sys.stderr, \
            "Matrix best matches column '%s' [%d:%d]." % (
            gs, nm, MATRIX.ncol())
        MAX_TO_SHOW = 5
        if len(missing) > MAX_TO_SHOW:
            print >>sys.stderr, \
                "Missing (showing %d of %d):" % (MAX_TO_SHOW, len(missing))
        else:
            print >>sys.stderr, "Missing from annotation:"
        for x_ in missing[:MAX_TO_SHOW]:
            print >>sys.stderr, x_
        raise AssertionError("I could not match the matrix to a geneset.")

    # Add the new column names to the MATRIX.
    MATRIX_new = MATRIX.matrix()
    name = arrayio.COL_ID
    if name not in MATRIX_new._col_names:
        name = MATRIX_new._synonyms[name]
    assert name in MATRIX_new._col_names, "I can not find the sample names."
    names = MATRIX_new.col_names(name)
    gs = genesets[0]
    genes = geneset2genes[gs]
    assert len(I_matrix) == len(I_geneset)
    for i in range(len(I_matrix)):
        names[I_matrix[i]] = "%s %s" % (
            names[I_matrix[i]], genes[I_geneset[i]])
    MATRIX_new._col_names[name] = names

    return MATRIX_new


def add_col_ids(MATRIX, geneset, ignore_missing):
    import os
    import sys
    import arrayio
    from genomicode import genesetlib
    from genomicode import matrixlib

    if not geneset:
        return MATRIX
    filename, genesets = _parse_file_gs(geneset)
    assert len(genesets) == 1

    # BUG: Should merge this code with append_col_ids.  Duplicated!
    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    ext = os.path.splitext(filename)[1].lower()
    for x in genesetlib.read_genesets(
            filename, allow_tdf=True, allow_duplicates=True):
        geneset, description, genes = x

        # Bug: sometimes will mis-identify TDF files as GMX.  The
        # first row will be interpreted as a description instead of a
        # gene (or annotation).  If the extension of the file isn't
        # gmx or gmt, and the description looks like a gene name, then
        # assume it's some sort of tdf file that has been
        # misidentified.
        #if not genesetlib._is_known_desc(description) and \
        #       ext not in [".gmx", ".gmt"]:
        if description and ext not in [".gmx", ".gmt"]:
            genes = [description] + genes

        geneset2genes[geneset] = genes
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Find an alignment between the sample names and the genesets.
    x = matrixlib.align_cols_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True)
    I_matrix, I_geneset, index = x
    if len(I_matrix) == 0:
        raise AssertionError("Matrixes doesn't match any gene sets.")
    elif len(I_matrix) != MATRIX.ncol() and not ignore_missing:
        # No matches.  Try to diagnose.
        gs = all_genesets[index]
        nm = len(I_matrix)
        missing = []
        col_names = MATRIX.col_names(arrayio.COL_ID)
        missing = [col_names[i] for i in range(len(col_names))
                   if i not in I_matrix]
        print >>sys.stderr, \
            "Matrix best matches column '%s' [%d:%d]." % (
            gs, nm, MATRIX.ncol())
        MAX_TO_SHOW = 5
        if len(missing) > MAX_TO_SHOW:
            print >>sys.stderr, \
                "Missing (showing %d of %d):" % (MAX_TO_SHOW, len(missing))
        else:
            print >>sys.stderr, "Missing from annotation:"
        for x_ in missing[:MAX_TO_SHOW]:
            print >>sys.stderr, x_
        raise AssertionError("I could not match the matrix to a geneset.")
    # END DUPLICATED SECTION.

    # Add the new column names to the MATRIX.
    MATRIX_new = MATRIX.matrix()

    gs = genesets[0]
    genes = geneset2genes[gs]
    assert len(I_matrix) == len(I_geneset)
    new_names = [None] * len(I_matrix)
    for i in range(len(I_matrix)):
        new_names[I_matrix[i]] = genes[I_geneset[i]]

    assert gs not in MATRIX_new._col_names
    MATRIX_new._col_names[gs] = new_names
    MATRIX_new._col_order.append(gs)
    return MATRIX_new


def reorder_col_indexes(MATRIX, indexes, count_headers):
    if not indexes:
        return MATRIX
    I_given = parse_indexes(MATRIX, False, indexes, count_headers)
    # Add back any indexes that weren't explicitly given by the user.
    I_other = [i for i in range(MATRIX.ncol()) if i not in I_given]
    I_all = I_given + I_other
    MATRIX_new = MATRIX.matrix(None, I_all)
    return MATRIX_new


def reorder_col_cluster(MATRIX, cluster, tree_file,
                        cluster_method, distance_method):
    from genomicode import cluster30
    from genomicode import clusterio
    from genomicode import matrixlib

    if not cluster:
        assert not tree_file
        return MATRIX
    if not MATRIX.nrow() or not MATRIX.ncol():
        return MATRIX

    cdata = cluster30.cluster_hierarchical(
        MATRIX, False, True, distance=distance_method, method=cluster_method)
    if tree_file:
        clusterio.write_atr_file(cdata.array_tree, open(tree_file, 'w'))

    # cluster30 will convert the MATRIX to PCL format, losing row
    # annotations.  Put missing row annotations back.
    all_annots = []
    for x in MATRIX.row_names():
        x = MATRIX.row_names(x)
        all_annots.append(x)
    x = matrixlib.align_rows_to_many_annots(
        cdata.matrix, all_annots, get_indexes=True)
    I_MATRIX, I_row_names, index = x
    assert I_MATRIX and I_row_names
    assert I_MATRIX == range(cdata.matrix.nrow())

    for header in MATRIX.row_names():
        if header in cdata.matrix.row_names():
            continue
        x = MATRIX.row_names(header)
        x = [x[i] for i in I_row_names]
        cdata.matrix._row_order.append(header)
        cdata.matrix._row_names[header] = x

    return cdata.matrix

    ## R = jmath.start_R()
    ## jmath.R_equals(MATRIX._X, "X")
    ## x = 'dist(t(X), method="%s")' % distance_method
    ## x = 'hclust(%s, method="%s")' % (x, cluster_method)
    ## I = list(R("%s$order" % x))   # 1-based indexes
    ## #I = list(R("hclust(dist(t(X)))$order"))   # 1-based indexes
    ## I = [x-1 for x in I]
    ## MATRIX_new = MATRIX.matrix(None, I)
    ## return MATRIX_new


def reorder_col_alphabetical(MATRIX, alphabetize):
    import arrayio
    from genomicode import jmath

    if not alphabetize:
        return MATRIX

    col_names = MATRIX.col_names(arrayio.COL_ID)
    I = jmath.order(col_names)
    MATRIX_new = MATRIX.matrix(None, I)
    return MATRIX_new


def reorder_col_byfile(MATRIX, filename, ignore_missing_cols):
    from genomicode import filelib

    if not filename:
        return MATRIX

    # File should be one line per sample name.
    samples = [x.strip() for x in filelib.openfh(filename)]
    I_row, I_col = MATRIX._index(col=samples)
    assert I_col, "I could not find any samples from %s." % filename
    if not ignore_missing_cols:
        assert len(I_col) == len(samples), "I could not find some samples."
    MATRIX_new = MATRIX.matrix(None, I_col)
    return MATRIX_new


def remove_duplicate_cols(MATRIX, filter_duplicate_cols):
    import arrayio

    if not filter_duplicate_cols:
        return MATRIX
    headers = MATRIX.col_names(arrayio.COL_ID)
    I = []
    seen = {}
    for i in range(len(headers)):
        if headers[i] in seen:
            continue
        seen[headers[i]] = 1
        I.append(i)
    x = MATRIX.matrix(None, I)
    return x


def remove_unnamed_cols(MATRIX, remove_cols):
    import arrayio

    if not remove_cols:
        return MATRIX
    headers = MATRIX.col_names(arrayio.COL_ID)
    I = []
    for i in range(len(headers)):
        if headers[i]:
            I.append(i)
    x = MATRIX.matrix(None, I)
    return x


def rename_duplicate_cols(MATRIX, rename_duplicate_cols):
    import arrayio

    if not rename_duplicate_cols:
        return MATRIX
    headers = MATRIX.col_names(arrayio.COL_ID)
    name2I = {}  # name -> list of indexes
    for i, name in enumerate(headers):
        if name not in name2I:
            name2I[name] = []
        name2I[name].append(i)

    nodup = headers[:]
    for (name, I) in name2I.iteritems():
        if len(I) < 2:
            continue
        for i in range(len(I)):
            nodup[I[i]] = "%s (%d)" % (name, i+1)

    MATRIX = MATRIX.matrix()  # don't change the original
    x = MATRIX._resolve_synonym(
        arrayio.COL_ID, MATRIX.col_names, MATRIX._synonyms)
    MATRIX._col_names[x] = nodup
    return MATRIX


## def remove_col_indexes(MATRIX, remove_col_indexes, count_headers):
##     if not remove_col_indexes:
##         return None
##     I = []
##     for indexes in remove_col_indexes:
##         x = parse_indexes(MATRIX, False, indexes, count_headers)
##         I.extend(x)

##     I_all = [x for x in range(MATRIX.ncol()) if x not in I]
##     return I_all


## def remove_col_ids(MATRIX, remove_col_ids):
##     import arrayio

##     if not remove_col_ids:
##         return MATRIX
##     x = ",".join(remove_col_ids)
##     x = x.split(",")
##     col_ids = [x.strip() for x in x]
##     names = MATRIX.col_names(arrayio.COL_ID)

##     I = []
##     not_found = {}.fromkeys(col_ids)
##     for i, name in enumerate(names):
##         if name in col_ids:
##             if name in not_found:
##                 del not_found[name]
##         else:
##             I.append(i)
##     assert not not_found, "I could not find: %s." % \
##         ", ".join(sorted(not_found))
##     x = MATRIX.matrix(None, I)
##     return x


def toupper_col_ids(MATRIX, toupper_col_ids):
    from arrayio import tab_delimited_format as tdf

    if not toupper_col_ids:
        return MATRIX
    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    MATRIX = MATRIX.matrix()
    x = MATRIX.col_names(tdf.SAMPLE_NAME)
    x = [x.upper() for x in x]
    MATRIX._col_names[tdf.SAMPLE_NAME] = x
    return MATRIX


def hash_col_ids(MATRIX, hash_col_ids):
    from arrayio import tab_delimited_format as tdf
    from genomicode import hashlib

    if not hash_col_ids:
        return MATRIX
    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    MATRIX = MATRIX.matrix()
    x = MATRIX.col_names(tdf.SAMPLE_NAME)
    x = [hashlib.hash_var(x) for x in x]
    MATRIX._col_names[tdf.SAMPLE_NAME] = x
    return MATRIX


def apply_re_col_ids(MATRIX, apply_re_col_ids):
    import re
    from arrayio import tab_delimited_format as tdf

    if not apply_re_col_ids:
        return MATRIX
    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    MATRIX = MATRIX.matrix()
    names = MATRIX.col_names(tdf.SAMPLE_NAME)
    for i in range(len(names)):
        m = re.search(apply_re_col_ids, names[i])
        if m:
            names[i] = m.group(1)
    MATRIX._col_names[tdf.SAMPLE_NAME] = names
    return MATRIX


def add_prefix_col_ids(MATRIX, prefix):
    from arrayio import tab_delimited_format as tdf

    if not prefix:
        return MATRIX
    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    MATRIX = MATRIX.matrix()
    names = MATRIX.col_names(tdf.SAMPLE_NAME)
    names = [prefix+x for x in names]
    MATRIX._col_names[tdf.SAMPLE_NAME] = names
    return MATRIX


def add_suffix_col_ids(MATRIX, suffix):
    from arrayio import tab_delimited_format as tdf

    if not suffix:
        return MATRIX
    if tdf.SAMPLE_NAME not in MATRIX.col_names():
        return MATRIX
    MATRIX = MATRIX.matrix()
    names = MATRIX.col_names(tdf.SAMPLE_NAME)
    names = [x+suffix for x in names]
    MATRIX._col_names[tdf.SAMPLE_NAME] = names
    return MATRIX


def _parse_tcga_barcode(barcode):
    # Return tuple of (patient, sample, aliquot, analyte).  sample,
    # aliquot, analyte may be None if these parts are not given in the
    # barcode.

    # Patient barcode  TCGA-02-0021
    # Sample barcode               -01B              + sample (without "-")
    # Aliquot barcode                  -02D          + portion
    # Analyte barcode                      -181-06   + plate and center
    assert len(barcode) >= 12, "Invalid barcode: %s" % barcode
    x = barcode.upper().split("-")
    assert len(x) >= 3, "Invalid barcode: %s" % barcode
    assert len(x) <= 7, "Invalid barcode: %s" % barcode

    sample = aliquot = analyte = None

    assert x[0] == "TCGA", "Invalid barcode: %s" % barcode
    assert len(x[1]) == 2, "Invalid barcode: %s" % barcode
    assert len(x[2]) == 4, "Invalid barcode: %s" % barcode
    patient = "-".join(x[:3])

    if len(x) >= 4:
        assert len(x[3]) >= 2 and len(x[3]) <= 3
        sample = x[3]
    if len(x) >= 5:
        assert len(x[4]) >= 2 and len(x[4]) <= 3
        aliquot = x[4]
    if len(x) >= 6:
        assert len(x) >= 7
        analyte = "-".join(x[6:8])

    return patient, sample, aliquot, analyte


def tcga_normal_only(MATRIX, cancer_only, ignore_non_tcga):
    from arrayio import tab_delimited_format as tdf

    if not cancer_only:
        return MATRIX
    assert tdf.SAMPLE_NAME in MATRIX.col_names()
    barcodes = MATRIX.col_names(tdf.SAMPLE_NAME)
    I = []
    for i, barcode in enumerate(barcodes):
        try:
            x = _parse_tcga_barcode(barcode)
        except AssertionError, x:
            # Keep all samples that don't look like a TCGA barcode.
            if ignore_non_tcga and str(x).startswith("Invalid barcode"):
                I.append(i)
                continue
            raise
        patient, sample, aliquot, analyte = x
        assert sample is not None, "sample missing from barcode"
        assert len(sample) >= 2
        sample = int(sample[:2])
        assert sample >= 1
        if sample >= 10 and sample < 20:
            I.append(i)
    x = MATRIX.matrix(None, I)
    return x


def tcga_primary_tumor_only(MATRIX, cancer_only, ignore_non_tcga):
    from arrayio import tab_delimited_format as tdf

    if not cancer_only:
        return MATRIX
    assert tdf.SAMPLE_NAME in MATRIX.col_names()
    barcodes = MATRIX.col_names(tdf.SAMPLE_NAME)
    I = []
    for i, barcode in enumerate(barcodes):
        try:
            x = _parse_tcga_barcode(barcode)
        except AssertionError, x:
            # Keep all samples that don't look like a TCGA barcode.
            if ignore_non_tcga and str(x).startswith("Invalid barcode"):
                I.append(i)
                continue
            raise
        patient, sample, aliquot, analyte = x
        assert sample is not None, "sample missing from barcode"
        assert len(sample) >= 2
        sample = int(sample[:2])
        assert sample >= 1
        if sample == 1:
            I.append(i)
    x = MATRIX.matrix(None, I)
    return x


def tcga_metastasis_only(MATRIX, cancer_only, ignore_non_tcga):
    from arrayio import tab_delimited_format as tdf

    if not cancer_only:
        return MATRIX
    assert tdf.SAMPLE_NAME in MATRIX.col_names()
    barcodes = MATRIX.col_names(tdf.SAMPLE_NAME)
    I = []
    for i, barcode in enumerate(barcodes):
        try:
            x = _parse_tcga_barcode(barcode)
        except AssertionError, x:
            # Keep all samples that don't look like a TCGA barcode.
            if ignore_non_tcga and str(x).startswith("Invalid barcode"):
                I.append(i)
                continue
            raise
        patient, sample, aliquot, analyte = x
        assert sample is not None, "sample missing from barcode"
        assert len(sample) >= 2
        sample = int(sample[:2])
        assert sample >= 1
        if sample == 6:  # also keep 7, additional metastatic?
            I.append(i)
    x = MATRIX.matrix(None, I)
    return x


def tcga_relabel_patient_barcodes(MATRIX, relabel, ignore_non_tcga):
    import arrayio
    if not relabel:
        return MATRIX

    name = arrayio.COL_ID
    if name not in MATRIX._col_names:
        name = MATRIX._synonyms[name]
    assert name in MATRIX._col_names, "I can not find the sample names."

    MATRIX_new = MATRIX.matrix()
    barcodes = MATRIX.col_names(name)
    for i in range(len(barcodes)):
        barcode = barcodes[i]
        try:
            x = _parse_tcga_barcode(barcode)
        except AssertionError, x:
            # Keep all samples that don't look like a TCGA barcode.
            if ignore_non_tcga and str(x).startswith("Invalid barcode"):
                pass
            else:
                raise
        else:
            barcode = x[0]
        barcodes[i] = barcode
    #x = [_parse_tcga_barcode(x) for x in barcodes]
    #x = [x[0] for x in x]
    MATRIX_new._col_names[name] = barcodes

    return MATRIX_new


def tcga_label_by_tissue_type(MATRIX, label_tissue, ignore_non_tcga):
    if not label_tissue:
        return MATRIX
    import arrayio

    name = arrayio.COL_ID
    if name not in MATRIX._col_names:
        name = MATRIX._synonyms[name]
    assert name in MATRIX._col_names, "I can not find the sample names."

    MATRIX_new = MATRIX.matrix()
    barcodes = MATRIX.col_names(name)
    
    for i in range(len(barcodes)):
        barcode = barcodes[i]
        try:
            x = _parse_tcga_barcode(barcode)
        except AssertionError, x:
            # Keep all samples that don't look like a TCGA barcode.
            if ignore_non_tcga and str(x).startswith("Invalid barcode"):
                pass
            else:
                raise
        else:
            assert len(x) >= 2
            patient, sample = x[:2]
            assert len(sample) >= 2
            sample_i = int(sample[:2])
            if sample_i == 1:
                barcode = "PRIMARY"
            elif sample_i == 2:
                barcode = "RECURRENT"
            elif sample_i == 6:
                barcode = "METASTATIC"
            elif sample_i == 7:
                barcode = "ADDITIONAL_METASTATIC"
            elif sample_i == 10:
                barcode = "NORMAL_BLOOD"
            elif sample_i == 11:
                barcode = "NORMAL_SOLID"
            else:
                raise AssertionError, "Unknown sample: %s" % sample
        barcodes[i] = barcode
    MATRIX_new._col_names[name] = barcodes

    return MATRIX_new


def select_row_indexes(MATRIX, indexes, count_headers):
    if not indexes:
        return None
    return parse_indexes(MATRIX, True, indexes, count_headers)


def select_row_ids(MATRIX, ids):
    # IDs is a list of gene IDs or names.  Combine them into a single string.
    ids = [x.strip() for x in ids]
    ids = ",".join(ids)
    if not ids:
        return None
    return parse_names(MATRIX, True, ids)


def select_row_string(MATRIX, row_annotation):
    # Format: <header>,<value>[,<value,...]
    if not row_annotation:
        return None

    x = _parse_file_string(row_annotation)
    header, values = x

    assert header in MATRIX.row_names(), "Missing header: %s" % header
    annots = MATRIX.row_names(header)
    I = []
    for i, annot in enumerate(annots):
        # Accepts matches for any of the values.
        if annot in values:
            I.append(i)
    return I


def select_row_numeric(MATRIX, row_annotation):
    # Format: <header>,<value>[,<value,...]
    if not row_annotation:
        return None

    x = _parse_file_num(row_annotation)
    header, values = x

    assert header in MATRIX.row_names(), "Missing header: %s" % header
    annots = MATRIX.row_names(header)
    I = []
    for i in range(len(annots)):
        x = float(annots[i])
        match = True
        for (modifier, value) in values:
            if modifier == "=":
                match = abs(x-value) < 1E-10
            elif modifier == "<":
                match = (x < value)
            elif modifier == ">":
                match = (x > value)
            elif modifier == "<=":
                match = (x <= value)
            elif modifier == ">=":
                match = (x >= value)
            else:
                raise AssertionError("Unknown modifier: %s" % modifier)
        #print x, modifier, value, match
        # Accepts matches for any of the values.
        if match:
            I.append(i)
    return I


def select_row_random(MATRIX, num_rows):
    import random
    if not num_rows:
        return None
    num_rows = int(num_rows)
    assert num_rows > 0 and num_rows < MATRIX.nrow()

    I = sorted(random.sample(range(MATRIX.nrow()), num_rows))
    return I


def select_row_genesets(MATRIX, genesets):
    if not genesets:
        return None
    I = []
    for geneset in genesets:
        x = parse_geneset(MATRIX, True, geneset)
        I.extend(x)
    # No duplicates.  Keep current order.
    i = 0
    while i < len(I):
        if I[i] in I[:i]:
            del I[i]
        else:
            i += 1
            
    return I


def select_row_annotation(MATRIX, row_annotation):
    # Format: <txt_file>,<header>,<value>[,<value,...]
    from genomicode import matrixlib

    if not row_annotation:
        return None

    x = _parse_file_annot(row_annotation)
    filename, header, values = x

    # Read the annotations.
    x = _read_annot_file(filename)
    header2annots, all_headers, all_annots = x

    # Align the annotations to the matrix file.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_annots, hash=True, get_indexes=True)
    I_matrix, I_annots, index = x
    assert len(I_matrix) == len(I_annots)
    #import arrayio
    #print index
    #for (m, a) in zip(I_matrix, I_annots):
    #    print m, a, MATRIX.row_names(arrayio.ROW_ID)[m], all_annots[index][a]

    annots = header2annots[header]
    I = []
    for (im, ia) in zip(I_matrix, I_annots):
        if annots[ia] in values:
            I.append(im)
    return I


def select_row_numeric_annotation(MATRIX, row_annotation):
    # Format: <txt_file>,<header>,<value>[,<value,...]
    from genomicode import matrixlib

    if not row_annotation:
        return None

    x = _parse_file_num_annot(row_annotation)
    filename, header, values = x

    # Read the annotations.
    x = _read_annot_file(filename)
    header2annots, all_headers, all_annots = x

    # Align the annotations to the matrix file.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_annots, hash=True, get_indexes=True)
    I_matrix, I_annots, index = x
    assert len(I_matrix) == len(I_annots)

    annots = header2annots[header]
    I = []
    for (im, ia) in zip(I_matrix, I_annots):
        x = float(annots[ia])
        match = True
        for (modifier, value) in values:
            if modifier == "=":
                match = (x == value)
            elif modifier == "<":
                match = (x < value)
            elif modifier == ">":
                match = (x > value)
            elif modifier == "<=":
                match = (x <= value)
            elif modifier == ">=":
                match = (x >= value)
            else:
                raise AssertionError("Unknown modifier: %s" % modifier)
        #print x, modifier, value, match
        # Accepts matches for any of the values.
        if match:
            I.append(im)
    return I


def select_row_nonempty(MATRIX, row_nonempty):
    # Format: <header>
    if not row_nonempty:
        return None

    header = row_nonempty
    assert header in MATRIX.row_names(), "Missing header: %s" % header
    annots = MATRIX.row_names(header)
    I = []
    for i, annot in enumerate(annots):
        if annot != "":
            I.append(i)
    return I


def select_row_maxvalue(MATRIX, maxvalue):
    if maxvalue is None:
        return None
    maxvalue = float(maxvalue)
    assert maxvalue >= 0 and maxvalue < 10000

    I = []  # indexes to keep.
    for i in range(len(MATRIX._X)):
        x = [x for x in MATRIX._X[i] if x is not None]
        if max(x) >= maxvalue:
            I.append(i)
    return I


def filter_row_minvalue(MATRIX, minvalue):
    if minvalue is None:
        return None
    minvalue = float(minvalue)
    assert minvalue >= 0 and minvalue < 50000

    I = []  # indexes to keep.
    for i in range(len(MATRIX._X)):
        x = [x for x in MATRIX._X[i] if x is not None]
        # Keep if all samples >= minvalue.
        if min(x) >= minvalue:
            I.append(i)
    return I


def filter_row_maxvalue(MATRIX, maxvalue):
    if maxvalue is None:
        return None
    maxvalue = float(maxvalue)
    assert maxvalue >= 0 and maxvalue < 50000

    I = []  # indexes to keep.
    for i in range(len(MATRIX._X)):
        x = [x for x in MATRIX._X[i] if x is not None]
        # Keep if all samples < maxvalue
        if max(x) < maxvalue:
            I.append(i)
    return I


def select_row_mean_value(MATRIX, mean_value):
    from genomicode import jmath
    
    if mean_value is None:
        return None
    mean_value = float(mean_value)
    assert mean_value >= 0 and mean_value < 10000

    means = jmath.mean(MATRIX._X)
    assert len(means) == len(MATRIX._X)
    I = []
    for i in range(len(means)):
        if means[i] >= mean_value:
            I.append(i)
    return I


def select_row_mean_var(MATRIX, filter_mean, filter_var):
    from genomicode import pcalib
    if filter_mean is None and filter_var is None:
        return None
    if filter_mean is not None:
        filter_mean = float(filter_mean)
    if filter_var is not None:
        filter_var = float(filter_var)

    assert filter_mean is None or (filter_mean >= 0 and filter_mean <= 1)
    assert filter_var is None or (filter_var >= 0 and filter_var <= 1)

    nrow = MATRIX.nrow()

    num_genes_mean = num_genes_var = None
    if filter_mean is not None:
        # Calculate the number of genes to keep.
        num_genes_mean = int(round((1.0 - filter_mean) * nrow))
    if filter_var is not None:
        # Calculate the number of genes to keep.
        num_genes_var = int(round((1.0 - filter_var) * nrow))
    I = pcalib.select_genes_mv(
        MATRIX._X, num_genes_mean, num_genes_var, test_for_missing_values=True)
    return I


def select_row_missing_values(MATRIX, perc_missing):
    # perc_missing of 0.25 means remove all rows with 25% or more
    # missing values.
    if perc_missing is None:
        return None
    assert perc_missing > 0 and perc_missing <= 1

    I = []  # indexes to keep.
    for i in range(len(MATRIX._X)):
        x = [x for x in MATRIX._X[i] if x is None]
        num_missing = len(x)
        pm = float(num_missing) / len(MATRIX._X[i])
        # If the percent missing is less than cutoff, then keep this
        # row.
        if pm < perc_missing:
            I.append(i)
    return I

def select_row_var(MATRIX, select_var):
    from genomicode import pcalib
    if select_var is None:
        return None
    select_var = int(select_var)
    assert select_var >= 1 and select_var <= MATRIX.nrow()
    I = pcalib.select_genes_var(MATRIX._X, select_var)
    #print select_var, len(I)
    return I


def select_row_delta(MATRIX, select_delta):
    if select_delta is None:
        return None
    select_delta = float(select_delta)
    assert select_delta > 0 and select_delta < 100000

    delta = []
    for x in MATRIX._X:
        delta.append(max(x) - min(x))

    I = [i for (i, x) in enumerate(delta) if x >= select_delta]
    return I


def select_row_fc(MATRIX, select_fc):
    import math
    if select_fc is None:
        return None
    select_fc = float(select_fc)
    assert select_fc > 0 and select_fc < 500

    log_fc = []
    for x in MATRIX._X:
        lfc = max(x) - min(x)
        log_fc.append(lfc)

    lfc_cutoff = math.log(select_fc, 2)
    I = [i for (i, x) in enumerate(log_fc) if x >= lfc_cutoff]
    return I


def select_row_num_samples_fc(MATRIX, num_samples, use_median=False):
    # By default, uses mean.  If use_median is true, then uses median.
    from genomicode import jmath

    if num_samples is None:
        return None
    num_samples = int(num_samples)
    assert num_samples > 0 and num_samples < MATRIX.ncol()

    # Calculate the means.
    if use_median:
        means = jmath.median(MATRIX._X)
    else:
        means = jmath.mean(MATRIX._X)
    assert len(means) == len(MATRIX._X)

    # In each row, count the number of samples that deviate at least
    # 2-fold (assuming log_2) from the mean.
    deviants = [None] * len(means)
    for i, row in enumerate(MATRIX._X):
        x = [abs(x-means[i]) for x in row]
        x = [x for x in x if x >= 1]
        deviants[i] = len(x)

    # Keep only the genes where at least num_samples deviate 2-fold
    # from the mean.
    I = [i for (i, x) in enumerate(deviants) if x >= num_samples]
    return I


def dedup_row_by_var(MATRIX, header):
    from genomicode import jmath
    if not header:
        return None

    assert header in MATRIX.row_names(), "Missing header: %s" % header

    annots = MATRIX.row_names(header)
    annot2i = {}  # annotation -> list of indexes
    for i, annot in enumerate(annots):
        if annot not in annot2i:
            annot2i[annot] = []
        annot2i[annot].append(i)

    assert_no_missing_values(MATRIX)
    variances = jmath.var(MATRIX._X)

    I = []
    for annot, indexes in annot2i.iteritems():
        if len(indexes) == 1:
            I.append(indexes[0])
            continue
        max_var = max_i = None
        for i in indexes:
            if max_var is None or variances[i] > max_var:
                max_var = variances[i]
                max_i = i
        assert max_i is not None
        I.append(max_i)
    I.sort()
    return I


def reverse_rows(MATRIX, reverse):
    if not reverse:
        return MATRIX
    I = range(MATRIX.nrow()-1, -1, -1)
    MATRIX_new = MATRIX.matrix(I, None)
    return MATRIX_new


def reorder_row_indexes(MATRIX, indexes, count_headers):
    if not indexes:
        return MATRIX
    I_given = parse_indexes(MATRIX, True, indexes, count_headers)
    # Add back any indexes that weren't explicitly given by the user.
    I_other = [i for i in range(MATRIX.nrow()) if i not in I_given]
    I_all = I_given + I_other
    MATRIX_new = MATRIX.matrix(I_all, None)
    return MATRIX_new


def reorder_row_cluster(
    MATRIX, cluster, tree_file, cluster_method, distance_method,
    indexes):
    from genomicode import cluster30
    from genomicode import clusterio

    if not cluster:
        assert not tree_file
        return MATRIX
    if not MATRIX.nrow() or not MATRIX.ncol():
        return MATRIX

    # Parse the indexes.
    if indexes:
        raise NotImplementedError

    cdata = cluster30.cluster_hierarchical(
        MATRIX, True, False, distance=distance_method,
        method=cluster_method)
    if tree_file:
        clusterio.write_gtr_file(cdata.gene_tree, open(tree_file, 'w'))
    return cdata.matrix

    ## R = jmath.start_R()
    ## jmath.R_equals(MATRIX._X, "X")
    ## x = 'dist(X, method="%s")' % distance_method
    ## x = 'hclust(%s, method="%s")' % (x, cluster_method)
    ## I = list(R("%s$order" % x))   # 1-based indexes
    ## I = [x-1 for x in I]
    ## MATRIX_new = MATRIX.matrix(I, None)
    ## return MATRIX_new


def reorder_row_cor(MATRIX, correlations, reverse_negative_cors,
                    indexes, count_headers):
    from genomicode import jmath

    if not correlations and not indexes:
        return MATRIX
    if not MATRIX.nrow() or not MATRIX.ncol():
        return MATRIX

    # Parse the correlation.
    vec = parse_correlations(correlations)

    # Parse the indexes.
    I = None
    if indexes:
        I = parse_indexes(MATRIX, False, indexes, count_headers)
    X = MATRIX.slice(None, I)

    #nrow, ncol = len(X), len(X[0])
    ncol = len(X[0])
    assert len(vec) == ncol, "vector %d indexes %d" % (len(vec), ncol)

    jmath.start_R()
    jmath.R_equals(X, "X")
    jmath.R_equals(vec, "vec")
    jmath.R("cors <- cor(vec, t(X))")
    if reverse_negative_cors:
        jmath.R("cors[cors < 0] <- -1 - cors[cors < 0]")
    O = list(jmath.R("order(cors, decreasing=TRUE)"))
    O = [x-1 for x in O]  # convert 1-based to 0-based indexes
    MATRIX_new = MATRIX.matrix(O, None)
    return MATRIX_new


def reorder_row_mean(MATRIX, means, indexes, count_headers):
    from genomicode import jmath

    if not means and not indexes:
        return MATRIX
    if not MATRIX.nrow() or not MATRIX.ncol():
        return MATRIX

    # Parse the indexes.
    I = None
    if indexes:
        I = parse_indexes(MATRIX, False, indexes, count_headers)
    X = MATRIX.slice(None, I)

    assert X, "empty matrix"
    nrow, ncol = len(X), len(X[0])
    assert nrow and ncol

    means = jmath.mean(X)
    O = jmath.order(means)
    MATRIX_new = MATRIX.matrix(O, None)
    return MATRIX_new


def rename_duplicate_rows(MATRIX, rename_duplicate_rows):
    import arrayio

    if not rename_duplicate_rows:
        return MATRIX
    names = MATRIX.row_names(arrayio.ROW_ID)
    name2I = {}  # name -> list of indexes
    for i, name in enumerate(names):
        if name not in name2I:
            name2I[name] = []
        name2I[name].append(i)

    nodup = names[:]
    for (name, I) in name2I.iteritems():
        if len(I) < 2:
            continue
        for i in range(len(I)):
            nodup[I[i]] = "%s_%d" % (name, i+1)

    MATRIX = MATRIX.matrix()  # don't change the original
    x = MATRIX._resolve_synonym(
        arrayio.ROW_ID, MATRIX.row_names, MATRIX._synonyms)
    MATRIX._row_names[x] = nodup
    return MATRIX


## def merge_rows_with_dup_values(MATRIX, merge_rows_with_dup_values):
##     if not merge_rows_with_dup_values:
##         return MATRIX

##     # While there are duplicate rows, merge them.
##     X = MATRIX._X
##     i = 0
##     while i < len(X):
##         dup = None
##         for j in range(i+1, len(X)):
##             delta = 0.0
##             for k in range(len(X[i])):
##                 delta += abs(X[i][k] - X[j][k])
##             if delta < 1E-5:
##                 dup = j
##                 break
##         if dup is None:
##             i += 1
##             continue
##         # Merge row i with row dup.
##         print "Dup", i, dup
##         print MATRIX.row_names("Ensembl Gene ID")[i]
##         print MATRIX.row_names("Ensembl Gene ID")[dup]
##         #_merge_rows(
##         raise NotImplementedError
##         pass
##     raise NotImplementedError



def align_rows(MATRIX, align_row_matrix, ignore_missing_rows):
    import os
    import arrayio

    if not align_row_matrix:
        return None
    assert os.path.exists(align_row_matrix), \
        "File not found: %s" % align_row_matrix

    ALIGN = arrayio.read(align_row_matrix)
    # Try all the headers and see if we can find a hit.
    # BUG: what if there's no header?
    best_I = []
    for header in ALIGN.row_names():
        ids = ALIGN.row_names(header)
        I_row, I_col = MATRIX._index(row=ids, row_header=arrayio.ROW_ID)
        # Bug: cannot just check length.  If there are duplicate rows,
        # the longest one may not be the one that matches the most
        # unique rows.
        if len(I_row) > len(best_I):
            best_I = I_row
        #I = I_row
        if len(best_I) == len(ids):
            break
    I = best_I
    if not ignore_missing_rows and len(ids) != len(I):
        # Diagnose problem here.
        x = ALIGN.row_names(arrayio.ROW_ID)
        ids_A = {}.fromkeys(x)
        x = MATRIX.row_names(arrayio.ROW_ID)
        ids_M = {}.fromkeys(x)
        missing = []  # In the align file, but not in my file.
        for id_ in ids_A:
            if id_ not in ids_M:
                missing.append(id_)
        if len(missing) < 10:
            for id_ in sorted(missing):
                print id_
        message = ("%d IDs from the align file are missing from the "
                   "matrix file." % len(missing))
        raise AssertionError(message)
    return I


def align_cols(MATRIX, align_col_matrix, ignore_missing_cols):
    import os
    import arrayio

    if not align_col_matrix:
        return None
    assert os.path.exists(align_col_matrix), \
        "File not found: %s" % align_col_matrix

    headers = MATRIX.col_names(arrayio.COL_ID)
    ALIGN = arrayio.read(align_col_matrix)
    # Try all the headers and see if we can find a hit.
    # Bug: what if there are duplicates in MATRIX or ALIGN?
    best_matches, best_I, best_header = None, [], ""
    for header in ALIGN.col_names():
        ids = ALIGN.col_names(header)
        I_row, I_col = MATRIX._index(col=ids, col_header=arrayio.COL_ID)
        # Count the number of unique matches.
        num_matches = {}
        for i in I_col:
            num_matches[headers[i]] = 1
        num_matches = len(num_matches)
        if best_matches is None or num_matches > best_matches:
            best_matches, best_I, best_header = num_matches, I_col, header
        if num_matches == len(ids):
            break
    I = best_I
    if not ignore_missing_cols and len(ids) != len(I):
        # Diagnose problem here.
        x = ALIGN.col_names(arrayio.COL_ID)
        ids_A = {}.fromkeys(x)
        x = MATRIX.col_names(arrayio.COL_ID)
        ids_M = {}.fromkeys(x)
        missing = []
        for id_ in ids_A:
            if id_ not in ids_M:
                missing.append(id_)
        if len(missing) < 10:
            for id_ in sorted(missing):
                print id_
        message = "I could not find %d column IDs." % len(missing)
        raise AssertionError(message)
    return I


def add_row_id(MATRIX, header):
    from genomicode import parselib

    if not header:
        return MATRIX

    assert header not in MATRIX.row_names(), "duplicate row header."

    MATRIX_new = MATRIX.matrix()
    x = ["GENE%s" % x for x in parselib.pretty_range(0, MATRIX.nrow())]
    MATRIX_new._row_order.insert(0, header)
    MATRIX_new._row_names[header] = x
    return MATRIX_new


def add_row_annot(MATRIX, row_annots, allow_unaligned_row_annot):
    # row_annot should be in the format <gmx/gmt_file>[,<geneset>].
    from genomicode import genesetlib
    from genomicode import matrixlib

    if not row_annots:
        return MATRIX
    filename, genesets = _parse_file_gs(row_annots)

    # Read all genesets out of the geneset file.
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    all_genes = []
    num_genes = None
    for x in genesetlib.read_genesets(
            filename, allow_tdf=True, allow_duplicates=True):
        geneset, description, genes = x
        geneset2genes[geneset] = genes
        if num_genes is None:
            num_genes = len(genes)
        assert len(genes) == num_genes, "%s %d %d" % (
            geneset, len(genes), num_genes)
        all_genesets.append(geneset)
        all_genes.append(genes)

    # Find an alignment between one of the matrix row_names and the
    # genesets.
    x = matrixlib.align_rows_to_many_annots(
        MATRIX, all_genes, hash=True, get_indexes=True, reorder_MATRIX=False)
    I_matrix, I_geneset, x = x
    assert not MATRIX.nrow() or len(I_matrix), \
        "I could not match the matrix to a geneset."

    MATRIX_new = MATRIX
    if not allow_unaligned_row_annot:
        MATRIX_new = MATRIX.matrix(I_matrix, None)

    # Add the new annotations to the MATRIX.
    for gs in genesets:
        assert gs in geneset2genes, "Missing geneset: %s" % gs
        assert gs not in MATRIX.row_names(), "duplicate name: %s" % gs

        genes = [""] * MATRIX.nrow()
        for (im, ig) in zip(I_matrix, I_geneset):
            genes[im] = geneset2genes[gs][ig]

        if not allow_unaligned_row_annot:
            # MATRIX_new can be subset of MATRIX.
            genes = [genes[i] for i in I_matrix]

        assert len(genes) == MATRIX_new.nrow()
        MATRIX_new._row_order.append(gs)
        MATRIX_new._row_names[gs] = genes

    return MATRIX_new


def remove_row_annot(MATRIX, name):
    assert name in MATRIX.row_names(), "I could not find name: %s" % name
    MATRIX_clean = MATRIX.matrix()
    Mc = MATRIX_clean
    assert name in Mc._row_names
    del Mc._row_names[name]
    if Mc._row_order and name in Mc._row_order:
        i = Mc._row_order.index(name)
        Mc._row_order.pop(i)
    # Bug: Does not remove from the synonyms list.
    return MATRIX_clean


def rename_row_annot(MATRIX, row_annot):
    if row_annot is None:
        return MATRIX
    assert "," in row_annot
    x = row_annot.split(",", 1)
    old_name, new_name = x
    #print MATRIX.row_names()
    assert old_name in MATRIX.row_names(), \
        "I could not find name: %s" % old_name
    assert new_name not in MATRIX.row_names(), \
        "Row name %s already exists." % new_name

    MATRIX_clean = MATRIX.matrix()
    Mc = MATRIX_clean
    assert old_name in Mc._row_names
    Mc._row_names[new_name] = Mc._row_names[old_name]
    del Mc._row_names[old_name]
    if Mc._row_order and old_name in Mc._row_order:
        i = Mc._row_order.index(old_name)
        Mc._row_order[i] = new_name
    if hasattr(Mc, "_synonyms"):
        synonyms = {}
        for (key, value) in Mc._synonyms.iteritems():
            if key == old_name:
                key = new_name
            if value == old_name:
                value = new_name
            synonyms[key] = value
        Mc.__dict__["_synonyms"] = synonyms
    return MATRIX_clean


def move_row_annot(MATRIX, move_row_annot):
    if move_row_annot is None:
        return MATRIX
    x = move_row_annot.split(",", 1)
    assert len(x) == 2
    old_index, new_index = x
    old_index = int(old_index)
    new_index = int(new_index)

    MATRIX_clean = MATRIX.matrix()
    assert MATRIX_clean._row_order
    order = MATRIX_clean._row_order[:]
    assert old_index >= 1 and old_index <= len(order), \
        "Invalid annotation index: %d" % old_index
    assert new_index >= 1 and new_index <= len(order), \
        "Invalid annotation index: %d" % new_index
    x = order.pop(old_index - 1)
    order.insert(new_index - 1, x)
    MATRIX_clean._row_order = order
    return MATRIX_clean


def concat_row_annot(MATRIX, concat_row_annot):
    if concat_row_annot is None:
        return MATRIX

    # concat_row_annot is in format:
    # <new_header>,<concat_char>,<header1>,<header2>,...
    x = concat_row_annot.split(",")
    assert len(x) >= 4
    new_header, concat_char = x[:2]
    headers = x[2:]

    assert new_header not in MATRIX.row_names()
    for header in headers:
        assert header in MATRIX.row_names(), "I could not find header: %s" % \
               header

    MATRIX_clean = MATRIX.matrix()
    assert MATRIX_clean._row_order

    old_annots = [MATRIX.row_names(x) for x in headers]
    new_annots = []
    for i in range(len(old_annots[0])):
        x = [x[i] for x in old_annots]
        x = concat_char.join(x)
        new_annots.append(x)

    MATRIX_clean._row_order.append(new_header)
    MATRIX_clean._row_names[new_header] = new_annots
    return MATRIX_clean


def _rlog_blind_h(X, outfile):
    from genomicode.jmath import start_R, R_fn, R_var, R_equals

    assert len(X), "empty matrix"
    assert len(X[0]), "empty matrix"
    ncol = len(X[0])
    assert ncol >= 2, "need at least 2 columns"
    
    R = start_R()
    R_fn("library", "DESeq2")
    R_equals(X, "X")
    # Make up some dummy design data.
    x = ["1"] + ["2"]*(ncol-1)
    R_fn("data.frame", condition=x, RETVAL="col.data")
    R_fn(
        "DESeqDataSetFromMatrix", countData=R_var("X"),
        colData=R_var("col.data"), design=R_var("~condition"), RETVAL="dds")
    R_fn("rlog", R_var("dds"), blind=R_var("TRUE"), RETVAL="x")
    R_fn("assay", R_var("x"), RETVAL="x")
    log_X_R = R["x"]

    # Convert this matrix into a Python object.  This matrix is
    # column-major.
    log_X_py = [[None]*log_X_R.ncol for i in range(log_X_R.nrow)]
    zzz = 0
    for j in range(log_X_R.ncol):
        for i in range(log_X_R.nrow):
            log_X_py[i][j] = log_X_R[zzz]
            zzz += 1

    # Write log_X_py to an outfile.
    handle = open(outfile, 'w')
    for x in log_X_py:
        print >>handle, "\t".join(map(str, x))
    handle.close()


def rlog_blind(X):
    # Fork a subprocess, because some R libraries generate garbage to
    # the screen.
    import os
    import tempfile

    outfile = None
    try:
        x, outfile = tempfile.mkstemp(dir="."); os.close(x)
        if os.path.exists(outfile):
            os.unlink(outfile)

        _run_forked(_rlog_blind_h, (X, outfile), {})

        assert os.path.exists(outfile), "failed"
        X_rlog = []
        for line in open(outfile):
            x = line.rstrip("\r\n").split("\t")
            x = map(float, x)
            X_rlog.append(x)
    finally:
        if FORKED_PID and outfile and os.path.exists(outfile):
            os.unlink(outfile)
    return X_rlog


def _calc_cpm_h(X, outfile):
    from genomicode.jmath import start_R, R_fn, R_var, R_equals

    assert len(X), "empty matrix"
    assert len(X[0]), "empty matrix"

    R = start_R()
    R_fn("library", "edgeR")
    R_equals(X, "X")
    R_fn("DGEList", counts=R_var("X"), RETVAL="dge")
    R_fn("cpm", R_var("dge"), RETVAL="x")
    log_X_R = R["x"]

    # Convert this matrix into a Python object.  This matrix is
    # column-major.
    log_X_py = [[None]*log_X_R.ncol for i in range(log_X_R.nrow)]
    zzz = 0
    for j in range(log_X_R.ncol):
        for i in range(log_X_R.nrow):
            log_X_py[i][j] = log_X_R[zzz]
            zzz += 1
            
    # Write log_X_py to an outfile.
    handle = open(outfile, 'w')
    for x in log_X_py:
        print >>handle, "\t".join(map(str, x))
    handle.close()
    

def calc_cpm(X):
    # Fork a subprocess, because some R libraries generate garbage to
    # the screen.
    import os
    import tempfile

    outfile = None
    try:
        x, outfile = tempfile.mkstemp(dir="."); os.close(x)
        if os.path.exists(outfile):
            os.unlink(outfile)

        _run_forked(_calc_cpm_h, (X, outfile), {})

        assert os.path.exists(outfile), "failed (%s)" % outfile
        X_rlog = []
        for line in open(outfile):
            x = line.rstrip("\r\n").split("\t")
            x = map(float, x)
            X_rlog.append(x)
    finally:
        # Only delete in parent.
        if FORKED_PID and outfile and os.path.exists(outfile):
            os.unlink(outfile)
    return X_rlog


def set_min_value(MATRIX, value):
    MATRIX = [x[:] for x in MATRIX]  # Make a copy.
    for i in range(len(MATRIX)):
        for j in range(len(MATRIX[i])):
            MATRIX[i][j] = max(MATRIX[i][j], value)
    return MATRIX


def normalize_rows_to(MATRIX, row_id):
    if row_id is None:
        return MATRIX
    I_row = parse_names(MATRIX, True, row_id)
    assert I_row, "Row not found: %s" % row_id
    assert len(I_row) == 1, "Row not unique: %s" % row_id
    I_norm = I_row[0]
    
    MATRIX = MATRIX.matrix()
    X = MATRIX._X
    x_norm = X[I_norm][:]
    for i in range(len(X)):
        for j in range(len(X[i])):
            X[i][j] = X[i][j] / x_norm[j]
    return MATRIX


def normalize_rows_sub(MATRIX, row_id):
    if row_id is None:
        return MATRIX
    I_row = parse_names(MATRIX, True, row_id)
    assert I_row, "Row not found: %s" % row_id
    assert len(I_row) == 1, "Row not unique: %s" % row_id
    I_norm = I_row[0]
    
    MATRIX = MATRIX.matrix()
    X = MATRIX._X
    x_norm = X[I_norm][:]
    for i in range(len(X)):
        for j in range(len(X[i])):
            X[i][j] = X[i][j] - x_norm[j]
    return MATRIX


def center_genes_mean(MATRIX, indexes):
    from genomicode import jmath

    I = []
    if indexes:
        I = parse_indexes(MATRIX, False, indexes, False)

    # Center the genes in place.
    assert_no_missing_values(MATRIX)
    X = MATRIX._X
    for i in range(len(X)):
        # Subtract the mean.
        X_i = X[i]
        X_sub = X_i
        if I:
            X_sub = [X_i[j] for j in I]
        m = jmath.mean(X_sub)
        X[i] = [x - m for x in X_i]


def center_genes_median(MATRIX, indexes):
    from genomicode import jmath

    I = []
    if indexes:
        I = parse_indexes(MATRIX, False, indexes, False)

    # Center the genes in place.
    assert_no_missing_values(MATRIX)
    X = MATRIX._X
    for i in range(len(X)):
        # Subtract the median.
        X_i = X[i]
        X_sub = X_i
        if I:
            X_sub = [X_i[j] for j in I]
        m = jmath.median(X_sub)
        X[i] = [x - m for x in X_i]


def normalize_genes_var(MATRIX, indexes):
    from genomicode import jmath

    I = []
    if indexes:
        I = parse_indexes(MATRIX, False, indexes, False)

    # Normalize the genes in place.
    assert_no_missing_values(MATRIX)
    X = MATRIX._X
    for i in range(len(X)):
        X_i = X[i]
        X_sub = X_i
        if I:
            X_sub = [X_i[j] for j in I]

        m = jmath.mean(X_sub)
        s = jmath.stddev(X_sub)

        # Subtract the mean.
        X_i = [x - m for x in X_i]
        # Normalize to stddev of 1.
        if s != 0:
            X_i = [x / s for x in X_i]
        # Add the mean back.
        X_i = [x + m for x in X_i]
        X[i] = X_i


def count_missing_values(MATRIX, header):
    if not header:
        return MATRIX

    # Count the number of missing values in each row.
    X = MATRIX._X
    missing_values = [None] * len(X)
    for i in range(len(X)):
        x = [x for x in X[i] if x is None]
        missing_values[i] = len(x)
    
    # Add the new annotations to the MATRIX.
    MATRIX_new = MATRIX.matrix()
    assert header not in MATRIX_new._row_names, "Duplicate header: %s" % header
    MATRIX_new._row_order.append(header)
    MATRIX_new._row_names[header] = missing_values

    return MATRIX_new


def median_fill_genes(MATRIX):
    from genomicode import jmath

    # Median-fill the genes in place.
    X = MATRIX._X
    for i in range(len(X)):
        m = jmath.safe_median(X[i])
        for j in range(len(X[i])):
            if X[i][j] is None:
                #print "FILLING", i, j
                X[i][j] = m


def zero_fill_genes(MATRIX):
    # Zero-fill the genes in place.
    X = MATRIX._X
    for i in range(len(X)):
        for j in range(len(X[i])):
            if X[i][j] is None:
                X[i][j] = 0.0


def impute_missing_values_knn(MATRIX, K):
    from genomicode import jmath
    if K is None:
        return MATRIX
    K = int(K)
    assert K >= 1 and K < 100

    #X = MATRIX._X
    jmath.start_R()
    jmath.R("library(impute)")
    jmath.R_equals(MATRIX._X, "X")
    jmath.R("obj <- impute.knn(X, k=%d)" % K)
    X_impute = jmath.R2py_matrix(jmath.R("obj$data"))
    MATRIX_new = MATRIX.matrix()
    MATRIX_new._X = X_impute
    return MATRIX_new
                

def add_missing_values(MATRIX, matrix_file):
    import os
    import arrayio
    
    if not matrix_file:
        return MATRIX
    assert os.path.exists(matrix_file)
    MATRIX_missing = arrayio.read(matrix_file)
    assert MATRIX.dim() == MATRIX_missing.dim()

    # Change MATRIX_new in place.
    MATRIX_new = MATRIX.matrix()
    X = MATRIX_new._X
    for i in range(len(X)):
        for j in range(len(X[i])):
            if MATRIX_missing._X[i][j] is None:
                X[i][j] = None
    return MATRIX_new
                

def _loess_normalize(X):
    from genomicode.jmath import start_R, R_fn, R_var, R_equals

    R = start_R()
    R_fn("library", "affy")
    R_equals(X, "X")
    R_fn(
        "normalize.loess", R_var("X"), RETVAL="X_norm",
        **{"log.it":R_var("FALSE")})
    X_norm_R = R["X_norm"]

    # Convert this matrix into a Python object.  This matrix is
    # column-major.
    X_norm_py = [[None]*X_norm_R.ncol for i in range(X_norm_R.nrow)]
    zzz = 0
    for j in range(X_norm_R.ncol):
        for i in range(X_norm_R.nrow):
            X_norm_py[i][j] = X_norm_R[zzz]
            zzz += 1
    return X_norm_py

def loess_normalize(MATRIX):
    import os
    import sys
    import tempfile

    outfile = pid = None
    try:
        x, outfile = tempfile.mkstemp(dir="."); os.close(x)
        if os.path.exists(outfile):
            os.unlink(outfile)

        # Fork a subprocess or R will generate garbage to the screen.
        r, w = os.pipe()
        pid = os.fork()

        if not pid:   # child
            os.close(r)
            w = os.fdopen(w, 'w')
            os.dup2(w.fileno(), sys.stdout.fileno())
            X_norm = _loess_normalize(MATRIX._X)
            handle = open(outfile, 'w')
            for x in X_norm:
                print >>handle, "\t".join(map(str, x))
            sys.exit(0)
        else:  # parent
            os.close(w)
            r = os.fdopen(r)
            r_output = r.read()  # output from R library
            os.waitpid(pid, 0)
            if not os.path.exists(outfile):
                assert os.path.exists(outfile), "%s\nR normalization failed" %\
                       r_output
            X_norm = []
            for line in open(outfile):
                cols = line.rstrip("\r\n").split("\t")
                cols = map(float, cols)
                X_norm.append(cols)
    finally:
        if pid:
            if outfile and os.path.exists(outfile):
                os.unlink(outfile)

    MATRIX_new = MATRIX.matrix()
    assert len(X_norm) == MATRIX.nrow()
    for x in X_norm:
        assert len(x) == MATRIX.ncol()
    MATRIX_new._X = X_norm
    return MATRIX_new


## def _match_rownames_to_geneset(MATRIX, all_genesets, geneset2genes):
##     # Return tuple of (I_matrix, I_geneset) or None if no match can be
##     # found.  Will find the largest match possible.

##     # Align every geneset to every row name in the matrix.
##     geneset_aligns = []  # list of (I_geneset, rowname, geneset)
##     matrix_aligns = []   # list of (I_matrix, rowname, geneset)
##     for name in MATRIX.row_names():
##         annots = MATRIX.row_names(name)
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             I_geneset = _align_geneset_to_matrix(annots, genes)
##             I_matrix = _align_matrix_to_geneset(annots, genes)
##             if I_geneset:
##                 x = I_geneset, name, gs
##                 geneset_aligns.append(x)
##             if I_matrix:
##                 x = I_matrix, name, gs
##                 matrix_aligns.append(x)

##     # First, try to find a geneset that matches the exactly matrix.
##     # Favor geneset_aligns over matrix_aligns to avoid changing the
##     # matrix.
##     for x in geneset_aligns:
##         I_geneset, rowname, geneset = x
##         I_matrix = range(MATRIX.nrow())
##         assert len(I_matrix) == len(I_geneset)
##         return I_matrix, I_geneset

##     # Otherwise, choose the match that generates the largest matrix.
##     I_matrix = None
##     for x in matrix_aligns:
##         I, rowname, geneset = x
##         if I_matrix is None or len(I) > len(I_matrix):
##             I_matrix = I
##     if I_matrix:
##         I_geneset = range(len(I_matrix))
##         return I_matrix, I_geneset

##     return None

## def _match_colnames_to_geneset(
##     MATRIX, all_genesets, geneset2genes, hash=False):
##     # Return I_geneset or None if no match can be found.  The sample
##     # names from MATRIX must match the names in the gene set exactly.
##     # MATRIX may be a subset of the gene set.
##     import arrayio
##     from genomicode import hashlib

##     annots = MATRIX.col_names(arrayio.COL_ID)
##     if hash:
##         annots = [hashlib.hash_R(x) for x in annots]
##         g2g = {}
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             genes = [hashlib.hash_R(x) for x in genes]
##             g2g[gs] = genes
##         geneset2genes = g2g

##     # Align every geneset to the col names in the matrix.
##     geneset_aligns = []  # list of (I_geneset, geneset)
##     for gs in all_genesets:
##         genes = geneset2genes[gs]
##         I_geneset = _align_geneset_to_matrix(annots, genes)
##         if I_geneset:
##             x = I_geneset, gs
##             geneset_aligns.append(x)

##     # Find a geneset that exactly matches the sample names in the matrix.
##     for x in geneset_aligns:
##         I_geneset, geneset = x
##         # Allow the geneset to be a superset of MATRIX.
##         if len(I_geneset) > MATRIX.ncol():
##             continue
##         return I_geneset

##     # Return None if not found.
##     return None

## def _best_match_colnames_to_geneset(
##     MATRIX, all_genesets, geneset2genes, hash=False):
##     # Return the (name, num_matches, list of matches, list of
##     # non-matches) for the geneset that best matches the colnames.
##     import arrayio
##     from genomicode import hashlib

##     annots = MATRIX.col_names(arrayio.COL_ID)
##     annots_orig = annots[:]
##     if hash:
##         annots = [hashlib.hash_R(x) for x in annots]
##         g2g = {}
##         for gs in all_genesets:
##             genes = geneset2genes[gs]
##             genes = [hashlib.hash_R(x) for x in genes]
##             g2g[gs] = genes
##         geneset2genes = g2g

##     # Align every geneset to the col names in the matrix.
##     best_gs = best_num_matches = best_matches = best_mismatches = None
##     for gs in all_genesets:
##         genes = geneset2genes[gs]
##         x = _num_matches_to_geneset(annots_orig, annots, genes)
##         num_matches, matches, mismatches = x
##         if best_num_matches is None or num_matches > best_num_matches:
##             best_gs = gs
##             best_num_matches = num_matches
##             best_matches = matches
##             best_mismatches = mismatches

##     return best_gs, best_num_matches, best_matches, best_mismatches

## def _num_matches_to_geneset(matrix_annots, hashed_annots, geneset_genes):
##     # Count the number of annotations that match the geneset.  Return
##     # number of matches, list of matches, list of mismatches.
##     assert len(matrix_annots) == len(hashed_annots)

##     count = 0
##     matches, mismatches = [], []
##     for i, annot in enumerate(hashed_annots):
##         if annot in geneset_genes:
##             count += 1
##             matches.append(matrix_annots[i])
##         else:
##             mismatches.append(matrix_annots[i])
##     return count, matches, mismatches

## def _align_geneset_to_matrix(matrix_annots, geneset_genes):
##     # Return a list of the indexes required to align the geneset to
##     # the rows of the matrix.  If it cannot be completely aligned,
##     # then return None.
##     from genomicode import jmath

##     if len(matrix_annots) > len(geneset_genes):
##         return None
##     I = jmath.match(matrix_annots, geneset_genes)
##     if None in I:
##         return None
##     return I

## def _align_matrix_to_geneset(matrix_annots, geneset_genes):
##     # Return a list of the indexes required to align the matrix to the
##     # geneset.  If it cannot be completely aligned, then return None.
##     from genomicode import jmath

##     if len(geneset_genes) > len(matrix_annots):
##         return None
##     I = jmath.match(geneset_genes, matrix_annots)
##     if None in I:
##         return None
##     return I

## def _find_col_header(MATRIX, col_names):
##     # Find the header in the MATRIX whose annotations can be found in
##     # col_names.  MATRIX can be a subset of col_names.  If multiple
##     # headers match, then just return the first one.  If none match,
##     # return None.
##     from genomicode import jmath
##     from genomicode import hashlib

##     h_col_names = [hashlib.hash_R(x) for x in col_names]

##     for name in MATRIX.col_names():
##         matrix_col_names = MATRIX.col_names(name)
##         h_matrix_col_names = [hashlib.hash_R(x) for x in matrix_col_names]
##         #I = jmath.match(h_col_names, h_matrix_col_names)
##         #if None in I:  # missing something
##         #    continue
##         I = jmath.match(h_matrix_col_names, h_col_names)
##         missing = []
##         for i in range(len(I)):
##             if I[i] is None:
##                 missing.append(matrix_col_names[i])
##         #print name, len(missing), missing
##         if len(missing):
##             continue
##         return name
##     return None


def _intersect_indexes(*indexes):
    # None means take all indexes.
    indexes = [x for x in indexes if x is not None]
    if not indexes:
        return None

    # Only want the indexes that occur in them all.  Preserve order.
    I = indexes[0]
    for x in indexes[1:]:
        I = [i for i in I if i in x]
    return I


## def _union_indexes(*indexes):
##     # Ignore None.  Does not preserve order.
##     indexes = [x for x in indexes if x is not None]

##     # Want indexes that occur in any of them.
##     I = []
##     for x in indexes:
##         I.extend(x)
##     I = sorted({}.fromkeys(I))
##     return I


def _dedup_indexes(I):
    # Get rid of duplicate indexes, preserving the original order of
    # the indexes.
    I = I[:]

    i = 0
    seen = {}
    while i < len(I):
        if I[i] in seen:
            I.pop(i)
        else:
            seen[I[i]] = 1
            i += 1
    return I


FORKED_PID = None
def _run_forked(fn, args, keywds):
    global FORKED_PID
    import os
    import sys

    r, w = os.pipe()
    pid = os.fork()
    FORKED_PID = pid
    if pid:   # Parent
        os.close(w)
        r = os.fdopen(r)
        output = r.read()
        #for line in r:
        #    pass
        #    #sys.stdout.write(line)   # output from R library
        os.waitpid(pid, 0)
    else:     # Child
        os.close(r)
        w = os.fdopen(w, 'w')
        os.dup2(w.fileno(), sys.stdout.fileno())
        fn(*args, **keywds)
        sys.exit(0)
    return output


def main():
    import sys
    import argparse
    import arrayio
    from genomicode import jmath
    from genomicode import matrixlib
    from genomicode import quantnorm

    parser = argparse.ArgumentParser(
        description="Slice the rows or columns of a matrix.")
    parser.add_argument("filename", nargs="+", help="Matrices to slice.")
    parser.add_argument(
        "--read_as_csv", default=False, action="store_true",
        help="Read as a CSV file.")
    parser.add_argument(
        "--skip_lines", default=None, type=int,
        help="Skip this number of lines in the file.")
    parser.add_argument(
        "--remove_comments",
        help="Remove rows that start with this character (e.g. '#')")
    parser.add_argument(
        "--num_header_cols", type=int,
        help="This number of columns are headers.  If not given, will guess.")
    parser.add_argument(
        "--clean_only", default=False, action="store_true",
        help="Only read_as_csv and remove_comments.")
    parser.add_argument(
        "--output_format", default="tdf", choices=["tdf", "gct"],
        help="Specify the format for the output file.")
    # If the user chooses an outfile, will need to implement it for
    # clean_only as well.
    #parser.add_argument(
    #    "-o", default=None, metavar="OUTFILE", dest="outfile",
    #    help="Save to this file.  By default, writes output to STDOUT.")

    group = parser.add_argument_group(
        title="Matrix Manipulation and Computations")
    group.add_argument(
        "--transpose",
        help="Transpose the matrix.  Format: <old row ID>,<new row ID>.  "
        "<old row ID> is the header of the column in the original file that "
        "should be used for the headers in the transposed file.  "
        "<new row ID> is what should be the name of the column of the IDs "
        "in the transposed file.")
    group.add_argument(
        "--transpose_nonmatrix", action="store_true",
        help="Just transpose the rows and columns.  May not be an "
        "expression matrix.  Should only have one file.  Ignores all "
        "other parameters.")
    group.add_argument(
        "--correlate", action="store_true",
        help="Calculate the pairwise correlation of the columns.")
    group.add_argument(
        "--correlate_against",
        help="The argument should be the name of another matrix file.  "
        "Correlate each row of this matrix to the row of the matrix "
        "in filename.  "
        "Each of the rows and columns should be aligned.")
    group.add_argument(
        "--correlate_some_vs_all", action="append",
        help="Correlate a selected list of rows against all other rows.  "
        "The argument should be a list of IDs (like --select_row_ids).  "
        "The output uses the first column as the row ID.")
    
    group.add_argument(
        "--calc_mean", action="store_true",
        help="Calculate the mean of each row.")
    group.add_argument(
        "--calc_sd", action="store_true",
        help="Calculate the standard deviation of each row.")
    group.add_argument(
        "--calc_range", action="store_true",
        help="Calculate the range (max - min) of each row.")
    group.add_argument(
        "--calc_pca_rows", action="store_true",
        help="Calculate the principal components of the rows.")

    group.add_argument(
        "--average_row_indexes", 
        help="Take the average of multiple rows of the matrix.  "
        "Add the mean as a new row in the matrix.  "
        "Format: 1-50,75 (1-based, inclusive).  "
        "Respects --row_indexes_include_headers.")
    
    group.add_argument(
        "--group_expression_by_samplename",
        help="Make a Prism formatted column table by grouping together "
        "the expression values.  All samples that share the same name "
        "will form a single group.  Format:<row ID to group>.")
    
    group = parser.add_argument_group(title="Missing values")
    group.add_argument(
        "--zerofill", action="store_true", help="Fill missing values with 0.")
    group.add_argument(
        "--fill_missing_values_median", action="store_true",
        help="Fill missing values with the median of the row.  Performed "
        "after logging, but before centering and normalizing.")
    group.add_argument(
        "--impute_missing_values_knn",
        help="Fill missing values using KNN.  Specify K, the number of "
        "neighbors to use (usually 10).")
    group.add_argument(
        "--count_missing_values",
        help="Count the number of missing values in each row and store in "
        "a new column.  The argument is the header for this column.")
    group.add_argument(
        "--add_missing_values",
        help="Put missing values back into a matrix.  The argument should "
        "be a filename that contains a matrix of the same dimensions with "
        "missing values.  I will add missing values to the same place.")

    group = parser.add_argument_group(title="Normalization")
    group.add_argument(
        "-l", "--log_transform", dest="log_transform", 
        action="store_true", help="Log transform the data.")
    group.add_argument(
        "--unlog", action="store_true",
        help="Exponentiate (unlog) transform the data.")
    group.add_argument(
        "--rlog_blind", action="store_true",
        help="Do a regularized log transformation (from DESeq2).  "
        "Use BLIND: most genes should not change expression across "
        "data set.")
    group.add_argument(
        "-q", "--quantile", "--qnorm", dest="quantile", action="store_true",
        help="Quantile normalize the data.")
    group.add_argument(
        "--cpm", action="store_true", help="Convert raw counts for NGS data "
        "into CPM (Uses edgeR package).")
    group.add_argument(
        "--loess", action="store_true", help="Loess normalize the data.")
    group.add_argument(
        "--gc", "--gene_center", dest="gene_center",
        choices=["mean", "median"], help="Center each gene by: mean, median.")
    group.add_argument(
        "--gc_subset_indexes",
        help="Will center the genes based on the mean (or median) of"
        "this subset of the samples.  Given as indexes, e.g. 1-5,8 "
        "(1-based, inclusive).")
    group.add_argument(
        "--gn", "--gene_normalize", dest="gene_normalize", default=None,
        choices=["ss", "var"],
        help="Normalize each gene by: ss (sum of squares), var (variance).")
    group.add_argument(
        "--gn_subset_indexes",
        help="Will normalize the genes based on the variance (or sum "
        "of squares) of this subset of the samples.  Given as indexes, "
        "e.g. 1-5,8 (1-based, inclusive).")
    group.add_argument(
        "--min_value", type=float,
        help="Set the minimum value for this matrix.  Done before logging.")
    group.add_argument(
        "--normalize_rows_to", 
        help="Normalize each row by dividing the values in a specified row.  "
        "The argument to this flag should be an ID that uniquely specifies "
        "the row to normalize to.")
    group.add_argument(
        "--normalize_rows_sub", 
        help="Normalize each row by subtracting a specified row.  "
        "The argument to this flag should be an ID that uniquely specifies "
        "the row to normalize to.")
    
    group = parser.add_argument_group(title="Column filtering")
    group.add_argument(
        "--select_col_indexes", default=[], action="append",
        help="Which columns to include e.g. 1-5,8 (1-based, inclusive).  "
        "(MULTI)")
    group.add_argument(
        "--col_indexes_include_headers", default=False, action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.  "
        "Applies to: select_col_indexes, reorder_col_indexes, "
        "reorder_row_cluster_subset_indexes, reorder_row_cor_subset_indexes")
    group.add_argument(
        "--select_col_ids", default=[], action="append",
        help="Comma-separate list of IDs to include.  (MULTI)")
    group.add_argument(
        "--select_col_annotation", default=[], action="append",
        help="Include only the cols where the annotation contains a "
        "specific value.  If multiple values are given, then will select "
        "the annotations that match any value (OR).  If this option is "
        "given multiple times, selects only the cols that match all the "
        "annotations (AND).  "
        "Format: <txt_file>,<header>,<value>[,<value,...].  (MULTI)")
    group.add_argument(
        "--select_col_genesets", default=[], action="append",
        help="Include only the samples from this geneset.  "
        "Format: <txt/gmx/gmt_file>,<geneset>[,<geneset>,...].  (MULTI)")
    group.add_argument(
        "--select_col_regex", default=[], action="append",
        help="Include columns that match this regular expression.  (MULTI)")
    group.add_argument(
        "--select_col_random", default=None, type=int,
        help="Select this number of columns at random.")
    group.add_argument(
        "--select_col_numeric_value",
        help="Include only the cols with a specific numeric value.  "
        "Format: <row_id>,<value>[,<value>,...].  "
        'If <value> starts with a "<", then will only find the rows where '
        "the annotation is less than <value>.  "
        'The analogous constraint will be applied for ">".  '
        "Accepts the match if any of the <value>s are true.")
    group.add_argument(
        "--reverse_col_selection", default=False, action="store_true",
        help="Remove the selected columns instead of keeping them.")
    ## group.add_argument(
    ##     "--remove_col_indexes", default=[], action="append",
    ##     help="Comma-separated list of indexes to remove.")
    ## group.add_argument(
    ##     "--remove_col_ids", default=[], action="append",
    ##     help="Comma-separated list of IDs to remove.")
    group.add_argument(
        "--remove_duplicate_cols", default=False, action="store_true",
        help="If a column is found multiple times, keep only the first one.")
    group.add_argument(
        "--remove_unnamed_cols", default=False, action="store_true",
        help="If a column has no name, remove it.")
    group.add_argument(
        "--reorder_col_indexes",
        help="Change the order of the data columns.  Give the indexes "
        "in the order that they should occur in the file, e.g. 1-5,8 "
        "(1-based, inclusive).  Can use --col_indexes_include_headers.")
    group.add_argument(
        "--reorder_col_cluster", default=False, action="store_true",
        help="Cluster the cols.")
    group.add_argument(
        "--col_tree_file", help="Write out the dendrogram of the col clusters "
        "in atr format.")
    # Cluster 3.0 implementation of clustering.
    group.add_argument(
        "--cluster_method", choices=[
            "complete", "single", "average", "centroid"],
        default="complete", help="Agglomeration method.  "
        "Applies to reorder_col_cluster and reorder_row_cluster.")
    group.add_argument(
        "--distance_method", choices=[
            "uncent-cor", "pearson", "abs-uncent-cor",
            "abs-pearson", "spearman", "kendall",
            "euclidean", "city-block"],
        default="euclidean", help="Distance measure for clustering.  "
        "Applies to reorder_col_cluster and reorder_row_cluster.")
    # For R implementation of clustering.
    #group.add_argument(
    #    "--cluster_method", choices=[
    #        "complete", "ward", "single", "average", "mcquitty",
    #        "median", "centroid"],
    #    default="complete", help="Agglomeration method.  "
    #    "Applies to reorder_col_cluster and reorder_row_cluster.")
    #group.add_argument(
    #    "--distance_method", choices=[
    #        "euclidean", "maximum", "manhattan", "canberra", "binary",
    #        "minkowski"],
    #    default="euclidean", help="Distance measure for clustering.  "
    #    "Applies to reorder_col_cluster and reorder_row_cluster.")
    group.add_argument(
        "--reorder_col_alphabetical", default=False, action="store_true",
        help="Sort the columns alphabetically.")
    group.add_argument(
        "--reorder_col_byfile",
        help="Reorder based on a file.  One line per sample name.")
    group.add_argument(
        "--align_col_matrix",
        help="Align the cols to a matrix in another file.")
    group.add_argument(
        "--ignore_missing_cols", default=False, action="store_true",
        help="Ignore any cols that can't be found in the align_col_matrix "
        "or reorder_col_byfile.")

    group = parser.add_argument_group(title="Column annotations")
    group.add_argument(
        "--toupper_col_ids", default=False, action="store_true",
        help="Convert column IDs to upper case.  "
        "(Done after relabel, remove, but before filtering duplicates.)")
    group.add_argument(
        "--hash_col_ids", default=False, action="store_true",
        help="Hash the column IDs to [a-zA-Z0-9_].")
    group.add_argument(
        "--rename_duplicate_cols", default=False, action="store_true",
        help="If multiple columns have the same header, make their names "
        "unique.")
    group.add_argument(
        "--rename_col_id", default=[], action="append",
        help="Rename a column ID.  Format: <from>,<to>.  "
        "<from> will be replaced with <to>.  (MULTI)")
    group.add_argument(
        "--replace_col_ids", default=[], action="append",
        help="Replace strings within the column IDs.  Format: <from>,<to>.  "
        "Instances of <from> will be replaced with <to>.  (MULTI)")
    group.add_argument(
        "--relabel_col_ids",
        help="Relabel the column IDs.  Format: <txt/gmx/gmt_file>,<geneset>.  "
        "One of the genesets in the file must match the current column IDs.")
    group.add_argument(
        "--append_col_ids",
        help="Append this to the end of the column IDs.  "
        "Format: <txt/gmx/gmt_file>,<geneset>.  "
        "One of the genesets in the file must match the current column IDs.")
    group.add_argument(
        "--add_col_ids", default=[], action="append",
        help="Add a new set of column IDs (MULTI).  "
        "Format: <txt/gmx/gmt_file>,<geneset>.  "
        "One of the genesets in the file must match the current column IDs.")
    group.add_argument(
        "--ignore_missing_labels", default=False, action="store_true",
        help="Any column labels that can't be found will not be relabeled.")
    group.add_argument(
        "--apply_re_col_ids",
        help="Apply a regular expression to the column IDs and take group 1.")
    group.add_argument(
        "--add_prefix_col_ids", help="Add a prefix to each column ID.")
    group.add_argument(
        "--add_suffix_col_ids", help="Add a suffix to each column ID.")

    group = parser.add_argument_group(title="TCGA barcode operations")
    group.add_argument(
        "--tcga_normal_only", default=False, action="store_true",
        help="Keep only the columns that contain normal sample.")
    group.add_argument(
        "--tcga_primary_tumor_only", default=False, action="store_true",
        help="Keep only the columns that contain primary solid tumor.")
    group.add_argument(
        "--tcga_metastasis_only", default=False, action="store_true",
        help="Keep only the columns that contain metastatic tumor.")
    group.add_argument(
        "--tcga_relabel_patient_barcodes", default=False, action="store_true",
        help="Sample names should be patient barcodes.")
    group.add_argument(
        "--tcga_label_by_tissue_type", default=False, action="store_true",
        help="Label PRIMARY, RECURRENT, METASTATIC, ADDITIONAL_METASTATIC, "
        "NORMAL_BLOOD, or NORMAL_SOLID.")
    group.add_argument(
        "--ignore_non_tcga", default=False, action="store_true",
        help="Keep all samples that don't look like a TCGA barcode.")

    group = parser.add_argument_group(title="Row filtering")
    group.add_argument(
        "--select_row_indexes",
        help="Which rows to include e.g. 1-50,75 (1-based, inclusive).")
    group.add_argument(
        "--row_indexes_include_headers", default=False, action="store_true",
        help="If not given (default), then row 1 is the first row "
        "with data.  If given, then row 1 is the very first row in "
        "the file, including the headers.")
    group.add_argument(
        "--select_row_ids", default=[], action="append",
        help="Comma-separated list of IDs (e.g. probes, gene names) "
        "to include.  (MULTI)")
    group.add_argument(
        "--select_row_string", default=[], action="append",
        help="Include only the rows where the columns contains a "
        "specific string value.  "
        "Format: <header>,<value>[,<value>,...].  "
        "Accepts the row if any of the <value>s match.  (MULTI)")
    group.add_argument(
        "--select_row_numeric", default=[], action="append",
        help="Include only the rows where the columns contains a "
        "specific numeric value.  "
        "Format: <header>,<value>[,<value>,...].  "
        'If <value> starts with a "<", then will only find the rows where '
        "the annotation is less than <value>.  "
        'The analogous constraint will be applied for ">".  '
        "Accepts the match if any of the <value>s are true.  (MULTI)")
    group.add_argument(
        "--select_row_random", help="Select this number of random rows.")
    group.add_argument(
        "--select_row_annotation",
        help="Include only the rows where the annotation contains a "
        "specific value.  Format: <txt_file>,<header>,<value>[,<value,...]")
    group.add_argument(
        "--select_row_numeric_annotation", default=[], action="append",
        help="Include only the rows where the annotation contains a "
        "specific numeric value.  "
        "Format: <txt_file>,<header>,<value>[,<value>,...].  "
        'If <value> starts with a "<", then will only find the rows where '
        "the annotation is less than <value>.  "
        'The analogous constraint will be applied for ">".  '
        "Accepts the match if any of the <value>s are true.  (MULTI)")
    group.add_argument(
        "--select_row_nonempty",
        help="Include only the rows that have a non-blank annotation.  "
        "Format: <header>")
    group.add_argument(
        "--select_row_maxvalue", default=None, type=float,
        help="Include the rows whose maximum value exceeds this.  "
        "E.g. if '5.0' is given, then keep the rows with at least one value "
        "greater than 5.")
    group.add_argument(
        "--filter_row_minvalue", default=None, type=float,
        help="Remove the rows where at least one sample has a value less "
        "than this.  E.g. if '5.0' is given, then remove the rows where "
        "any sample is less than 5.")
    group.add_argument(
        "--filter_row_maxvalue", default=None, type=float,
        help="Remove the rows where at least one sample has a value more "
        "than this.  E.g. if '500.0' is given, then remove the rows where "
        "any sample is more than 500.")
    group.add_argument(
        "--select_row_genesets", default=[], action="append",
        help="Include only the IDs from this geneset.  "
        "Format: <txt/gmx/gmt_file>,<geneset>[,<geneset>,...].  (MULTI)")
    group.add_argument(
        "--filter_row_by_mean", default=None, type=float,
        help="Remove this percentage of rows that have the lowest mean.  "
        "Should be between 0 and 1.")
    group.add_argument(
        "--filter_row_by_var", default=None, type=float,
        help="Remove this percentage of rows that have the lowest variance.  "
        "Should be between 0 and 1.")
    group.add_argument(
        "--filter_row_by_missing_values", default=None, type=float,
        help="Remove the rows that has at least this percent of missing "
        "bvalues.  e.g. 0.25 means remove all rows with 25%% or more missing "
        "values.")
    ## group.add_argument(
    ##     "--merge_rows_with_dup_values", default=False, action="store_true",
    ##     help="Merge the annotations of the rows whose values are duplicated.")
    group.add_argument(
        "--dedup_row_by_var",
        help="If multiple rows have the same annotation, select the one "
        "with the highest variance.  The value of this parameter should "
        "be the header of the column that contains duplicate annotations.")
    group.add_argument(
        "--select_row_mean_value", default=None, type=float,
        help="Keep only the rows where the mean is at least this number.")
    group.add_argument(
        "--select_row_var", default=None, type=int,
        help="Keep this number of rows with the highest variance.")
    group.add_argument(
        "--select_row_delta", default=None, type=float,
        help="Keep only the rows where the difference between the maximum "
        "and minimum value is at least this.")
    group.add_argument(
        "--select_row_fc", default=None, type=float,
        help="Keep only the rows with at least this fold change between "
        "highest and lowest sample (assuming log_2 values).")
    group.add_argument(
        "--select_row_num_samples_fc_mean", default=None, type=int,
        help="Keep only the rows where at least this number of samples "
        "deviate at least 2 fold change from the mean "
        "(assuming log_2 values).")
    group.add_argument(
        "--select_row_num_samples_fc_median", default=None, type=int,
        help="Keep only the rows where at least this number of samples "
        "deviate at least 2 fold change from the median "
        "(assuming log_2 values).")
    group.add_argument(
        "--reverse_rows", default=False, action="store_true",
        help="Reverse the order of the rows.")
    group.add_argument(
        "--reorder_row_indexes",
        help="Change the order of the data rows.  Give the indexes "
        "in the order that they should occur in the file, e.g. 1-5,8 "
        "(1-based, inclusive).  Can use --row_indexes_include_headers.")
    group.add_argument(
        "--reorder_row_cluster", default=False, action="store_true",
        help="Cluster the rows.")
    group.add_argument(
        "--reorder_row_cluster_subset_indexes",
        help="Will cluster the rows based on a subset of the samples.")

    group.add_argument(
        "--row_tree_file", help="Write out the dendrogram of the row clusters "
        "in gtr format.")
    group.add_argument(
        "--reorder_row_cor",
        help="Reorder the rows based on a correlation to this vector.  "
        "This should be a comma-separated list of numbers, e.g. "
        "0,0,0,1,1,1")
    group.add_argument(
        "--reverse_negative_cors", default=False, action="store_true",
        help="UNDOCUMENTED")
    group.add_argument(
        "--reorder_row_cor_subset_indexes",
        help="Will reorder rows based on correlation to this subset of "
        "samples.  Can use --col_indexes_include_headers.")
    group.add_argument(
        "--reorder_row_mean", action="store_true",
        help="Reorder the rows based on mean expression.")
    group.add_argument(
        "--reorder_row_mean_subset_indexes", 
        help="Reorder the rows based on the mean expression of this "
        "subset of samples.  Can use --col_indexes_include_headers.")
    group.add_argument(
        "--align_row_matrix",
        help="Align the rows to a matrix in another file.")
    group.add_argument(
        "--ignore_missing_rows", default=False, action="store_true",
        help="Ignore any rows that can't be found in the align_row_matrix.")

    group = parser.add_argument_group(title="Row annotations")
    group.add_argument(
        "--add_row_id",
        help="Add a unique row ID.  This should be the name of the header.")
    group.add_argument(
        "--add_row_annot", action="append", default=[],
        help="Add a geneset as a new annotation for the matrix.  "
        "The format should be: <txt/gmx/gmt_file>,<geneset>[,<geneset>].  "
        "Each geneset in the file should contain the same number of "
        "genes as the matrix.  One of the genesets should be align-able "
        "to the IDs of this matrix.  (MULTI)")
    group.add_argument(
        "--allow_unaligned_row_annot", default=False, action="store_true",
        help="If the matrix contains rows not in the annotation file, "
        "fill them with empty annotations (rather than dropping the "
        "row).")
    group.add_argument(
        "--remove_row_annot", action="append", default=[],
        help="Remove this annotations from the matrix.  (MULTI)")
    group.add_argument(
        "--rename_row_annot", action="append", default=[],
        help="Rename this header.  "
        "The format should be: <old_name>,<new_name>.  (MULTI)")
    group.add_argument(
        "--move_row_annot", action="append", default=[],
        help="Move this header.  "
        "The format should be: <old_index>,<new_index>.  "
        "The indexes are 1-based.  (MULTI)")
    group.add_argument(
        "--concat_row_annot",
        help="Concatenate multiple row annotations.  "
        "The format should be: <new_header>,<concat_char>,<header>,...")
    group.add_argument(
        "--rename_duplicate_rows", default=False, action="store_true",
        help="If multiple rows have the same ID, make their names "
        "unique.")

    args = parser.parse_args()
    assert len(args.filename) >= 1

    if args.num_header_cols is not None:
        assert args.num_header_cols > 0 and args.num_header_cols < 100

    if args.transpose_nonmatrix:
        assert len(args.filename) == 1
        transpose_nonmatrix(args.filename[0])
        return

    x = read_matrices(
        args.filename, args.skip_lines, args.read_as_csv, args.remove_comments,
        args.clean_only, args.num_header_cols)
    fmt_module, matrices = x
    if len(matrices) == 1:
        MATRIX = matrices[0]
    else:
        # Merge the matrices into one big file.
        matrices = matrixlib.align_rows(*matrices)
        MATRIX = matrixlib.merge_matrices(*matrices)
    if not MATRIX.nrow():
        return

    MATRIX = transpose_matrix(MATRIX, args.transpose)
    MATRIX = correlate_matrix(MATRIX, args.correlate)
    MATRIX = correlate_against_matrix(MATRIX, args.correlate_against)
    MATRIX = correlate_some_vs_all(MATRIX, args.correlate_some_vs_all)
    MATRIX = calc_mean(MATRIX, args.calc_mean)
    MATRIX = calc_sd(MATRIX, args.calc_sd)
    MATRIX = calc_range(MATRIX, args.calc_range)
    MATRIX = calc_pca_rows(MATRIX, args.calc_pca_rows)
    MATRIX = average_row_indexes(
        MATRIX, args.average_row_indexes, args.row_indexes_include_headers)
    MATRIX = group_expression_by_samplename(
        MATRIX, args.group_expression_by_samplename)


    # Slice to a submatrix.
    I01 = select_row_indexes(
        MATRIX, args.select_row_indexes, args.row_indexes_include_headers)
    I02 = select_row_ids(MATRIX, args.select_row_ids)
    x = [select_row_string(MATRIX, annot) for annot in args.select_row_string]
    I03 = _intersect_indexes(*x)
    x = [select_row_numeric(MATRIX, annot)
         for annot in args.select_row_numeric]
    I04 = _intersect_indexes(*x)
    I05 = select_row_random(MATRIX, args.select_row_random)
    I06 = select_row_genesets(MATRIX, args.select_row_genesets)
    I07 = select_row_annotation(MATRIX, args.select_row_annotation)
    x = [select_row_numeric_annotation(MATRIX, annot)
         for annot in args.select_row_numeric_annotation]
    I08 = _intersect_indexes(*x)
    I09 = select_row_nonempty(MATRIX, args.select_row_nonempty)
    I10 = select_row_maxvalue(MATRIX, args.select_row_maxvalue)
    I11 = filter_row_minvalue(MATRIX, args.filter_row_minvalue)
    I12 = filter_row_maxvalue(MATRIX, args.filter_row_maxvalue)
    I13 = select_row_mean_value(MATRIX, args.select_row_mean_value)
    I14 = select_row_mean_var(
        MATRIX, args.filter_row_by_mean, args.filter_row_by_var)
    I15 = select_row_var(MATRIX, args.select_row_var)
    I16 = select_row_delta(MATRIX, args.select_row_delta)
    I17 = select_row_fc(MATRIX, args.select_row_fc)
    I18 = select_row_num_samples_fc(
        MATRIX, args.select_row_num_samples_fc_mean, use_median=False)
    I19 = select_row_num_samples_fc(
        MATRIX, args.select_row_num_samples_fc_median, use_median=True)
    I20 = select_row_missing_values(MATRIX, args.filter_row_by_missing_values)
    I_row = _intersect_indexes(
        I01, I02, I03, I04, I05, I06, I07, I08, I09, I10, I11, I12, I13,
        I14, I15, I16, I17, I18, I19, I20)

    I01 = select_col_indexes(
        MATRIX, args.select_col_indexes, args.col_indexes_include_headers)
    #I2 = remove_col_indexes(
    #    MATRIX, args.remove_col_indexes, args.col_indexes_include_headers)
    I03 = select_col_ids(MATRIX, args.select_col_ids)
    I04 = select_col_genesets(MATRIX, args.select_col_genesets)
    I05 = select_col_annotation(MATRIX, args.select_col_annotation)
    I06 = []
    for regex in args.select_col_regex:
        I06.extend(select_col_regex(MATRIX, regex))
    I06 = sorted({}.fromkeys(I06))
    if not args.select_col_regex:
        I06 = None   # if not given, use all columns.
    I07 = select_col_random(MATRIX, args.select_col_random)
    I08 = select_col_numeric_value(MATRIX, args.select_col_numeric_value)
    I_col = _intersect_indexes(I01, I03, I04, I05, I06, I07, I08)
    if args.reverse_col_selection:
        I_col = [i for i in range(MATRIX.ncol()) if i not in I_col]
    MATRIX = MATRIX.matrix(I_row, I_col)

    # Reorder the rows and columns by indexes.  Do this before
    # removing columns.  Do this before adding or removing
    # annotations.
    MATRIX = reorder_row_indexes(
        MATRIX, args.reorder_row_indexes, args.row_indexes_include_headers)
    MATRIX = reorder_col_indexes(
        MATRIX, args.reorder_col_indexes, args.col_indexes_include_headers)
    MATRIX = reorder_col_alphabetical(MATRIX, args.reorder_col_alphabetical)

    ## # Merge the rows with duplicated values.
    ## MATRIX = merge_rows_with_dup_values(
    ##     MATRIX, args.merge_rows_with_dup_values)

    # Remove row annotations.
    for name in args.remove_row_annot:
        MATRIX = remove_row_annot(MATRIX, name)

    # Add a unique row ID.
    MATRIX = add_row_id(MATRIX, args.add_row_id)

    # Add row annotations.
    for annot in args.add_row_annot:
        MATRIX = add_row_annot(MATRIX, annot, args.allow_unaligned_row_annot)

    # Rename the row annotations.  Do this after removing and adding.
    for x in args.rename_row_annot:
        MATRIX = rename_row_annot(MATRIX, x)

    # Move the row annotations.  Do this at the end.
    for x in args.move_row_annot:
        MATRIX = move_row_annot(MATRIX, x)

    MATRIX = concat_row_annot(MATRIX, args.concat_row_annot)

    # Relabel the column IDs.
    MATRIX = rename_col_id(
        MATRIX, args.rename_col_id, args.ignore_missing_labels)
    MATRIX = replace_col_ids(MATRIX, args.replace_col_ids)
    MATRIX = relabel_col_ids(
        MATRIX, args.relabel_col_ids, args.ignore_missing_labels)
    MATRIX = append_col_ids(
        MATRIX, args.append_col_ids, args.ignore_missing_labels)
    for col_ids in args.add_col_ids:
        MATRIX = add_col_ids(MATRIX, col_ids, args.ignore_missing_labels)

    # Remove col IDs.  Do this after relabeling.
    #MATRIX = remove_col_ids(MATRIX, args.remove_col_ids)

    # Convert col IDs to upper case.  Do after relabeling and
    # removing, but before filtering duplicates.
    MATRIX = toupper_col_ids(MATRIX, args.toupper_col_ids)
    MATRIX = apply_re_col_ids(MATRIX, args.apply_re_col_ids)
    MATRIX = add_prefix_col_ids(MATRIX, args.add_prefix_col_ids)
    MATRIX = add_suffix_col_ids(MATRIX, args.add_suffix_col_ids)
    MATRIX = hash_col_ids(MATRIX, args.hash_col_ids)

    # Filter TCGA columns.
    MATRIX = tcga_normal_only(
        MATRIX, args.tcga_normal_only, args.ignore_non_tcga)
    MATRIX = tcga_primary_tumor_only(
        MATRIX, args.tcga_primary_tumor_only, args.ignore_non_tcga)
    MATRIX = tcga_metastasis_only(
        MATRIX, args.tcga_metastasis_only, args.ignore_non_tcga)
    MATRIX = tcga_relabel_patient_barcodes(
        MATRIX, args.tcga_relabel_patient_barcodes, args.ignore_non_tcga)
    MATRIX = tcga_label_by_tissue_type(
        MATRIX, args.tcga_label_by_tissue_type, args.ignore_non_tcga)

    # Filter after relabeling.
    MATRIX = remove_duplicate_cols(MATRIX, args.remove_duplicate_cols)
    MATRIX = remove_unnamed_cols(MATRIX, args.remove_unnamed_cols)

    # Rename duplicate rows and columns.
    MATRIX = rename_duplicate_rows(MATRIX, args.rename_duplicate_rows)
    MATRIX = rename_duplicate_cols(MATRIX, args.rename_duplicate_cols)

    # Reorder based on a file.
    MATRIX = reorder_col_byfile(
        MATRIX, args.reorder_col_byfile, args.ignore_missing_cols)

    # Align to the align_file.  Do this as close to the end as
    # possible, after everything else removed and added.
    I_row = align_rows(MATRIX, args.align_row_matrix, args.ignore_missing_rows)
    MATRIX = MATRIX.matrix(I_row, None)
    I_col = align_cols(MATRIX, args.align_col_matrix, args.ignore_missing_cols)
    MATRIX = MATRIX.matrix(None, I_col)

    if args.min_value is not None:
        MATRIX._X = set_min_value(MATRIX._X, args.min_value)
    MATRIX = normalize_rows_to(MATRIX, args.normalize_rows_to)
    MATRIX = normalize_rows_sub(MATRIX, args.normalize_rows_sub)
    
    # Log transform, if requested.  Do before other changes to
    # expression values (quantile, center, normalize).
    if args.log_transform:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)
    if args.unlog:
        MATRIX._X = jmath.exp(MATRIX._X, base=2)
    if args.rlog_blind:
        MATRIX._X = rlog_blind(MATRIX._X)

    if args.cpm:
        MATRIX._X = calc_cpm(MATRIX._X)

    # Median fill.  Do after logging, but before quantile, centering,
    # and normalizing.
    MATRIX = count_missing_values(MATRIX, args.count_missing_values)
    if args.zerofill:
        zero_fill_genes(MATRIX)
    if args.fill_missing_values_median:
        median_fill_genes(MATRIX)
    MATRIX = impute_missing_values_knn(MATRIX, args.impute_missing_values_knn)
    MATRIX = add_missing_values(MATRIX, args.add_missing_values)
        

    # Quantile normalize, if requested.  After log.
    if args.quantile:
        MATRIX = quantnorm.normalize(MATRIX)

    # Loess normalize, if requested.  After log.
    if args.loess:
        MATRIX = loess_normalize(MATRIX)

    # This has to happen before any normalization by variance, but
    # after logging and quantile normalization.
    I = dedup_row_by_var(MATRIX, args.dedup_row_by_var)
    MATRIX = MATRIX.matrix(I, None)

    # Preprocess the expression values.
    if args.gene_center == "mean":
        center_genes_mean(MATRIX, args.gc_subset_indexes)
    elif args.gene_center == "median":
        center_genes_median(MATRIX, args.gc_subset_indexes)

    if args.gene_normalize == "ss":
        raise NotImplementedError
    elif args.gene_normalize == "var":
        normalize_genes_var(MATRIX, args.gn_subset_indexes)

    # Cluster the rows and columns.  Do this after normalizing, zero-fill, log.
    MATRIX = reorder_row_cluster(
        MATRIX, args.reorder_row_cluster, args.row_tree_file,
        args.cluster_method, args.distance_method,
        args.reorder_row_cor_subset_indexes)
    MATRIX = reorder_col_cluster(
        MATRIX, args.reorder_col_cluster, args.col_tree_file,
        args.cluster_method, args.distance_method)

    # Reorder the rows based on correlation.
    MATRIX = reorder_row_cor(
        MATRIX, args.reorder_row_cor, args.reverse_negative_cors,
        args.reorder_row_cor_subset_indexes, args.col_indexes_include_headers)

    # Reorder the rows based on correlation.
    MATRIX = reorder_row_mean(
        MATRIX, args.reorder_row_mean, 
        args.reorder_row_mean_subset_indexes, args.col_indexes_include_headers)

    # Reverse the rows.  Do after all the selection.  Do after
    # aligning to a file.
    MATRIX = reverse_rows(MATRIX, args.reverse_rows)

    # Write the outfile (in the same format).
    handle = sys.stdout
    #if args.outfile:
    #    handle = open(args.outfile, 'w')

    # Cannot always write in the same format.  For example, if you add
    # annotations that aren't handled by that format.  To be safe,
    # convert to a TDF (by default) and write that out.
    assert args.output_format in ["tdf", "gct"]
    to_format = arrayio.tdf
    if args.output_format == "gct":
        to_format = arrayio.gct_format
    MATRIX = arrayio.convert(MATRIX, to_format=to_format)
    to_format.write(MATRIX, handle)


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
