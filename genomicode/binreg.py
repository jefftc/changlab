"""

Functions:
binreg
binreg_raw

find_binreg_20
format_data_files
format_pref_file
format_exec_file
format_predictions

is_logged_array_data
log_matrix_if_needed

align_rows             # WAS align_matrices
align_cols
are_rows_aligned       # is_matrices_aligned
are_cols_aligned
describe_unaligned_rows

read_matrices
merge_gct_matrices
strip_affx_control_probes

Classes:
BinregParams

"""
# Private functions:
# _hash_geneid
# _hash_sampleid

import os, sys


class BinregParams:
    def __init__(
        self, binreg_version=None, cross_validate=None, make_plots=None,
        num_genes=None, num_metagenes=None,
        quantile_normalize=None, shift_scale_normalize=None,
        num_burnin=None, num_iterations=None, num_skips=None,
        credible_interval=None):
        # Set defaults
        if binreg_version is None:
            binreg_version = 2
        if cross_validate is None:
            cross_validate = 1
        if make_plots is None:
            make_plots = 0
        if num_genes is None:
            num_genes = 100
        if num_metagenes is None:
            num_metagenes = 2
        if quantile_normalize is None:
            quantile_normalize = 1
        if shift_scale_normalize is None:
            shift_scale_normalize = 1
        if num_burnin is None:
            num_burnin = 1000
        if num_iterations is None:
            num_iterations = 5000
        if num_skips is None:
            num_skips = 1
        if credible_interval is None:
            credible_interval = 95

        # Make sure inputs are valid.
        assert binreg_version in [1, 2]
        assert cross_validate in [0, 1]
        assert make_plots in [0, 1]
        assert num_genes > 0
        assert num_metagenes > 0
        assert quantile_normalize in [0, 1]
        assert shift_scale_normalize in [0, 1]
        assert num_burnin >= 0
        assert num_iterations > 0
        assert num_skips >= 0
        assert credible_interval >= 0 and credible_interval <= 100

        self.binreg_version = binreg_version
        self.cross_validate = cross_validate
        self.make_plots = make_plots
        self.num_genes = num_genes
        self.num_metagenes = num_metagenes
        self.quantile_normalize = quantile_normalize
        self.shift_scale_normalize = shift_scale_normalize
        self.num_burnin = num_burnin
        self.num_iterations = num_iterations
        self.num_skips = num_skips
        self.credible_interval = credible_interval

def binreg(train0, train1, test, is_logged, params, 
           matlab=None, binreg_path=None, outpath=None):
    outpath = outpath or "."

    desc_file = os.path.join(outpath, "description.txt")
    express_file = os.path.join(outpath, "expression.txt")
    x = format_data_files(train0, train1, test)
    open(desc_file, 'w').write(x[0])
    open(express_file, 'w').write(x[1])

    x = binreg_raw(
        express_file, desc_file, is_logged, params, matlab=matlab,
        binreg_path=binreg_path, outpath=outpath)
    return x

def binreg_raw(expression_file, description_file, is_logged, params,
               matlab=None, binreg_path=None, outpath=None):
    # expression_file contains a gene x sample matrix of the
    # expression values.  There should be no headers.  The first row
    # contains 0/1/2 depending on whether the sample is train0,
    # train1, or test.  The description_file contains a list of the
    # gene names or probeset IDs and should be parallel to the
    # expression_file.
    import subprocess
    
    # Set defaults.
    matlab = matlab or "matlab"
    outpath = outpath or "."
    binreg_path = find_binreg_20(binreg_path)
    pref_file = os.path.join(outpath, "preferences.dat")

    assert binreg_path is not None, "I could not find Binreg 2.0"
    binreg_path = os.path.realpath(binreg_path)
    pref_file = os.path.realpath(pref_file)

    # Check parameters.
    assert os.path.exists(expression_file), "Missing: %s" % expression_file
    assert os.path.exists(description_file), "Missing: %s" % description_file
    assert is_logged in [0, 1]
    assert binreg_path, "I could not find a complete BinReg 2.0 distribution."
    
    x = format_pref_file(expression_file, description_file, is_logged, params)
    open(pref_file, 'w').write(x)

    # Run Binreg from Matlab.
    cwd = os.getcwd()
    try:
        os.chdir(outpath)
        matlab_args = [
            "-nosplash", "-nodesktop", "-nodisplay", "-nojvm"]
        x = " ".join(matlab_args)
        cmd = "%s %s" % (matlab, x)
        #print cmd
        #w, r = os.popen4(cmd, bufsize=0)
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout

        x = format_exec_file(binreg_path, pref_file)
        w.write(x)
        w.close()
    finally:
        os.chdir(cwd)
    return r

def find_binreg_20(default_path):
    # Return the path to BinReg2.0 or None.
    search_paths = [
        default_path, 
        "BinReg2.0",
        "/Volumes/Users/jchang/projects/geclassify/data/BinReg2.0",
        "/home/jchang/projects/geclassify/data/BinReg2.0",
        "/home/jefftc/projects/geclassify/data/BinReg2.0",
        "/home/vis/jefftc/projects/geclassify/data/BinReg2.0",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        # Test for some well-known Binreg files.
        files = [
            "binreg.m", "binreg_batch.m", "Mbinregsvd.m", "README.txt",
            "quantnorm.m", "std_rows_to_y.m", 
            "Binreg2.0.Matlab.tutorial.ppt"]
        complete = True   # Is this distribution complete.
        for file in files:
            filename = os.path.join(spath, file)
            if not os.path.exists(filename):
                complete = False
                break
        if not complete:
            continue
        path = spath
        break
    return path

def format_data_files(train0, train1, test):
    # test can be None.
    from StringIO import StringIO
    import arrayio

    assert train0.nrow() > 0 and train0.ncol() > 0
    assert train1.nrow() > 0 and train1.ncol() > 0
    assert not test or (test.nrow() > 0 and test.ncol() > 0)
    assert are_rows_aligned(train0, train1), "matrices not aligned"
    assert not test or are_rows_aligned(train0, test), "matrices not aligned"

    X_train0 = train0.value()
    X_train1 = train1.value()
    if test:
        X_test = test.value()

    # Merge the matrices.
    X_all = []
    x = [0]*len(X_train0[0]) + [1]*len(X_train1[0])
    if test:
        x = x + [2]*len(X_test[0])
    X_all.append(x)
    for i in range(len(X_train0)):
        x = X_train0[i] + X_train1[i]
        if test:
            x = x + X_test[i]
        X_all.append(x)

    # Write the description file.
    ids = train0.row_names(arrayio.ROW_ID)
    desc_handle = StringIO()
    for id in ids:
        print >>desc_handle, id
    desc_handle.seek(0)
    desc_str = desc_handle.read()

    # Write the expression_file
    express_handle = StringIO()
    for x in X_all:
        print >>express_handle, "\t".join(map(str, x))
    express_handle.seek(0)
    express_str = express_handle.read()

    return desc_str, express_str

def format_pref_file(expression_file, description_file, is_logged, params):
    # Return a string with the formatted preference file.
    expression_file = os.path.realpath(expression_file)
    description_file = os.path.realpath(description_file)
    x = [
        expression_file,
        description_file,
        params.num_genes,
        params.num_metagenes,
        params.num_burnin,
        params.num_iterations,
        params.num_skips,
        params.credible_interval,
        is_logged,
        params.quantile_normalize,
        params.shift_scale_normalize,
        params.binreg_version,
        params.cross_validate,
        params.make_plots,
        ]
    x = map(str, x)
    x = "\n".join(x) + "\n"
    return x

def format_exec_file(binreg_path, pref_file):
    from StringIO import StringIO

    handle = StringIO()
    w = handle.write
    
    # The BinReg gamrnd needs to mask the Matlab builtin, so add the
    # binreg path to the beginning.  The functions give different
    # results.
    w("addpath '%s' -begin;\n" % binreg_path)
    #w("addpath '%s' -end;\n" % binreg_path)
    w("binreg_batch('%s');\n" % pref_file)
    w("quit;\n")
    handle.seek(0)
    return handle.read()

def format_predictions(train0, train1, test, outpath=None):
    # Return a string.  test can be None.
    from StringIO import StringIO
    import arrayio
    import filelib

    outpath = outpath or "."

    opj = os.path.join
    fitted_file = opj(outpath, "trainingcases.txt")
    xval_file = opj(outpath, "crossvalidation.txt")
    predict_file = opj(outpath, "validationcases.txt")
    
    assert os.path.exists(fitted_file), "Missing trainingcases.txt."
    assert os.path.exists(xval_file), "Missing crossvalidation.txt."
    if test:
        assert os.path.exists(predict_file), "Missing validationcases.txt."

    samples = train0.col_names(arrayio.COL_ID)+train1.col_names(arrayio.COL_ID)
    if test:
        samples = samples + test.col_names(arrayio.COL_ID)
        
    format = "index:d type:d prob:f lower_ci:f upper_ci:f mgene:f"
    d_fit = [d for d in filelib.read_row(fitted_file, format)]
    d_xval = [d for d in filelib.read_row(xval_file, format)]
    d_pred = []
    if test:
        d_pred = [d for d in filelib.read_row(predict_file, format)]
    assert len(d_fit) == len(d_xval)
    # crossvalidation.txt + validationcases.txt
    assert len(d_xval)+len(d_pred) == len(samples)

    # Make the types pretty.
    type_str = ["train0", "train1", "test"]
    for d in d_fit+d_xval+d_pred:
        assert d.type in [0, 1, 2], "Unknown type: %d" % d.Type
        d.type = type_str[d.type]

    # Add the method.
    for d in d_xval+d_pred:
        d.method = "PREDICTED"
    for d in d_fit:
        d.method = "FITTED"

    # Add the sample names.
    for i in range(len(d_fit)):
        d_fit[i].sample = samples[i]
        d_xval[i].sample = samples[i]
    for i in range(len(d_pred)):
        d_pred[i].sample = samples[i+len(d_fit)]

    # Indexes here are 0-based.
    handle = StringIO()
    header = ["Index", "Sample", "Type", "Method", "Probability",
              "Lower CI", "Upper CI", "Metagene"]
    print >>handle, "\t".join(header)
    
    for i, d in enumerate(d_fit+d_xval+d_pred):
        x = (d.index-1, d.sample, d.type, d.method,
             d.prob, d.lower_ci, d.upper_ci, d.mgene)
        print >>handle, "\t".join(map(str, x))
    handle.seek(0)
    return handle.read()

def is_logged_array_data(data):
    CUTOFF = 1000

    # Count the number of scores greater than the cutoff.
    nrow, ncol = data.dim()
    X = data.slice()
    total = nrow * ncol
    num_greater = 0
    for i in range(nrow):
        for j in range(ncol):
            if X[i][j] >= CUTOFF:
                num_greater += 1

    # If all scores < CUTOFF, then is logged.
    if num_greater == 0:
        return True
    # If a portion of the scores >= CUTOFF, then is not logged.
    #if num_greater >= 10 or float(num_greater)/total >= 0.10:
    if num_greater >= 10 or float(num_greater)/total >= 0.01:
        return False
    # If only a few scores >= CUTOFF, then maybe there's an error in the file.
    raise AssertionError, "Detected some outliers [%d:%d].  Log?" % (
        num_greater, total)

def log_matrix_if_needed(DATA):
    import jmath
    if is_logged_array_data(DATA):
        return DATA
    DATA = DATA.matrix()
    DATA._X = jmath.log(DATA._X, base=2, safe=1)
    return DATA

def align_rows(*matrices):
    """Aligns matrices by ROW_ID.  Return a list of the matrices after
    the rows are aligned.  Raises an exception if no rows are common
    to all matrices.

    """
    import arrayio

    header = arrayio.ROW_ID
    if not matrices:
        return []
    for m in matrices:
        assert m.nrow() > 0 and m.ncol() > 0

    # Get the intersection of the IDs.
    all_ids = [None]*len(matrices) # matrix index -> list of hids
    for i in range(len(matrices)):
        all_ids[i] = _hash_many_geneids(matrices[i].row_names(header))
    
    ids_hash = {}  # hid -> matrix_i -> row_i, id
    for i, ids in enumerate(all_ids):
        for j, id in enumerate(ids):
            hid = _hash_geneid(id)
            # If I have not seen this hashed ID before, then add it to
            # the dictionary.
            if hid not in ids_hash:
                ids_hash[hid] = {}
            # If a dataset has duplicate IDs, then take the first one only.
            if i in ids_hash[hid]:
                continue
            ids_hash[hid][i] = j, id

    # Use only the IDs that occur in all data files.  Use the order of
    # the samples from first data set.
    ids = [id for id in all_ids[0] if len(ids_hash[id]) == len(matrices)]
    assert len(ids) > 0, "The data sets all have different row IDs."
            
    # Align the rows by the ids.
    aligned = [None] * len(matrices)
    for i in range(len(matrices)):
        I = [ids_hash[id][i][0] for id in ids]
        x = matrices[i].matrix(row=I)
        aligned[i] = x

    assert are_rows_aligned(*aligned), "matrices not aligned"

    return aligned

def align_cols(*matrices):
    """Aligns matrices by COL_ID.  Return a list of the matrices after
    the columns are aligned.  Raises an exception if no columns are
    common to all matrices.

    """
    import arrayio

    header = arrayio.COL_ID
    if not matrices:
        return []
    for m in matrices:
        assert m.nrow() > 0 and m.ncol() > 0

    # Get the intersection of the IDs.
    all_ids = [None]*len(matrices) # matrix index -> list of hids
    for i in range(len(matrices)):
        all_ids[i] = _hash_many_sampleids(matrices[i].col_names(header))
    
    ids_hash = {}  # hid -> matrix_i -> col_i, id
    for i, ids in enumerate(all_ids):
        for j, id in enumerate(ids):
            hid = _hash_sampleid(id)
            # If I have not seen this hashed ID before, then add it to
            # the dictionary.
            if hid not in ids_hash:
                ids_hash[hid] = {}
            # If a dataset has duplicate IDs, then take the first one only.
            if i in ids_hash[hid]:
                continue
            ids_hash[hid][i] = j, id
            
    # Use only the IDs that occur in all data files.  Use the order of
    # the samples from first data set.
    ids = [id for id in all_ids[0] if len(ids_hash[id]) == len(matrices)]
    assert len(ids) > 0, "The data sets have different column IDs."

    # Align the columns by the ids.
    aligned = [None] * len(matrices)
    for i in range(len(matrices)):
        I = [ids_hash[id][i][0] for id in ids]
        x = matrices[i].matrix(col=I)
        aligned[i] = x

    assert are_cols_aligned(*aligned), "matrices not aligned"

    return aligned

def are_rows_aligned(*matrices):
    import arrayio

    header = arrayio.ROW_ID
    
    if len(matrices) <= 1:
        return True

    # Check the number of rows in each matrix.
    for i in range(1, len(matrices)):
        if matrices[0].nrow() != matrices[i].nrow():
            return False

    # Check the names of the rows.
    hnames = [_hash_many_geneids(x.row_names(header)) for x in matrices]
    for i in range(1, len(matrices)):
        if hnames[0] != hnames[i]:
            return False
    return True

def describe_unaligned_rows(*matrices):
    # Return a text string that describes where the rows are not aligned.
    import arrayio
    import parselib

    header = arrayio.ROW_ID
    
    if len(matrices) <= 1:
        return "Only 1 matrix.  Must be aligned."

    # Check the number of rows in each matrix.
    num_rows = [x.nrow() for x in matrices]
    if min(num_rows) != max(num_rows):
        x = "Matrices have differing number of rows: %s" % ", ".join(
            map(parselib.pretty_int, num_rows))
        return x

    # Check the names of the rows.
    names = [x.row_names(header) for x in matrices]
    hnames = [_hash_many_geneids(x.row_names(header)) for x in matrices]
    bad_rows = []
    for i in range(matrices[0].nrow()):
        unaligned =False
        for j in range(1, len(matrices)):
            if hnames[0][i] != hnames[j][i]:
                unaligned = True
        if unaligned:
            x = [names[j][i] for j in range(len(matrices))]
            x = "Row %s: %s" % (parselib.pretty_int(i+1), ", ".join(x))
            bad_rows.append(x)
    if not bad_rows:
        return "Matrices are aligned."

    total_bad = len(bad_rows)
    if total_bad > 10:
        bad_rows = bad_rows[:10]
        bad_rows.append("...")
    x = "%s of %s rows are unaligned." % (
        parselib.pretty_int(total_bad),
        parselib.pretty_int(matrices[0].nrow()))
    lines = [x] + bad_rows
    return "\n".join(lines)

def are_cols_aligned(*matrices):
    import arrayio

    header = arrayio.COL_ID

    if len(matrices) <= 1:
        return True

    # Check the number of columns in each matrix.
    for i in range(1, len(matrices)):
        if matrices[0].ncol() != matrices[i].ncol():
            return False

    # Check the names of the columns.
    hnames = [_hash_many_sampleids(x.col_names(header)) for x in matrices]
    for i in range(1, len(matrices)):
        if hnames[0] != hnames[i]:
            return False
    return True

def read_matrices(filenames):
    """Read a list of matrices and align them.  filenames is a list of
    the matrix files to read.  Returns a tuple where the first element
    is a list of the matrices read, and the second is the aligned
    matrix.

    """
    import arrayio
    import filelib

    for filename in filenames:
        assert filelib.exists(filename)

    # Load the files.
    DATA = []
    for filename in filenames:
        try:
            x = arrayio.read(filename)
        except (SystemError, KeyboardInterrupt, MemoryError), x:
            raise
        except Exception, x:
            # Can diagnose which file failed here.
            # raise
            raise Exception, "Problem reading %s: %s" % (
                repr(filename), str(x))
        DATA.append(x)

    for d, filename in zip(DATA, filenames):
        f = os.path.split(filename)[1]
        #print "%s has %d genes and %d samples." % (f, d.nrow(), d.ncol())

    # Align the matrices.
    ALIGNED = align_rows(*DATA)
    #if DATA:
    #    print "Merged file has %d genes." % DATA[0].nrow()
    #sys.stdout.flush()

    return DATA, ALIGNED

def merge_gct_matrices(*matrices):
    import arrayio
    import Matrix

    assert len(matrices)
    assert are_rows_aligned(*matrices)
    for m in matrices:
        assert arrayio.gct_format.is_matrix(m)

    # Get the values in list of lists format.
    X = [x.value() for x in matrices]
    X_all = []
    for i in range(len(X[0])):   # each row
        x = []
        for j in range(len(X)):  # each matrix
            x.extend(X[j][i])
        X_all.append(x)

    # The sample names may not be unique.
    row_order = ["NAME", "DESCRIPTION"]
    col_order = []
    row_names = {}
    col_names = {}
    synonyms = {}

    # Should be in GCT format.
    matrix0 = matrices[0]
    assert len(matrix0.row_names()) == 2, "Invalid file format."
    # Just assume the first column is the NAME and second is the
    # DESCRIPTION.  Allow more flexibility in the actual column
    # headers.
    name, description = matrix0.row_names()
    row_names["NAME"] = matrix0.row_names(name)
    row_names["DESCRIPTION"] = matrix0.row_names(description)

    for m in matrices:
        assert len(m.col_names()) == 1
    col_order = matrix0.col_names()
    assert len(col_order) == 1

    x = []
    for m in matrices:
        x.extend(m.col_names(m.col_names()[0]))
    col_names[col_order[0]] = x

    synonyms[arrayio.ROW_ID] = "NAME"
    synonyms[arrayio.COL_ID] = col_order[0]

    DATA = Matrix.InMemoryMatrix(
        X_all, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    DATA = Matrix.add_synonyms(DATA, synonyms)
    assert arrayio.gct_format.is_matrix(DATA)
    return DATA

def strip_affx_control_probes(DATA, psid_header=None):
    # psid_header should be the header that contains the Probe Set IDs.
    import arrayio

    header = psid_header or arrayio.ROW_ID
    
    ID = _hash_many_geneids(DATA.row_names(header))
    ID_at = [x for x in ID if x.endswith("_at")]
    ID_good = [x for x in ID if not x.startswith("affx")]
    assert ID, "No IDs."
    assert ID_at, "IDs do not look like Affymetrix IDs."
    assert ID_good, "I could not find any non-control probes."

    # Get the non-AFFX IDs.
    DATA_s = DATA.matrix(row=ID_good, row_header=header)

    return DATA_s

def _hash_geneid(id):
    return id.strip().lower()

def _hash_many_geneids(ids):
    return [_hash_geneid(x) for x in ids]

def _hash_sampleid(id):
    # Hash the sample names so that small differences are ignored.  R
    # will change the sample names, so if one data set has been
    # through R but not the other, the names will be different.
    # 2554_6933_32492_Mock1_HG-U133A+2
    # X2554_6933_32492_Mock1_HG.U133A.2
    import re

    x = id

    # If there are alphanumeric characters, then assume that
    # punctuation isn't meaningful.
    if re.search(r"\w", x):
        # Change all non-words to '.', like R does.  (This does not
        # change underscores.)
        x = re.subn(r"\W", ".", x)[0]

    # Ignore initial X.
    if re.match(r"X[\d\w]", x):
        x = x[1:]

    # Make case insensitive.
    x = x.lower()

    return x

def _hash_many_sampleids(ids):
    return [_hash_sampleid(x) for x in ids]
