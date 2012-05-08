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

## Moved to matrixlib.
#align_rows             # WAS align_matrices
#align_cols
#are_rows_aligned       # is_matrices_aligned
#are_cols_aligned
#describe_unaligned_rows
#
#read_matrices
#merge_gct_matrices
#merge_matrices

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
    import config
    
    search_paths = [
        default_path, 
        config.binreg20_path,
        "BinReg2.0",
        "/Volumes/Users/jchang/projects/geclassify/data/BinReg2.0",
        "/home/jchang/projects/geclassify/data/BinReg2.0",
        "/home/jefftc/projects/geclassify/data/BinReg2.0",
        "/home/vis/jefftc/projects/geclassify/data/BinReg2.0",
        "/home/changlab/SIGNATURE/BinReg2.0",
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
    import matrixlib

    assert train0.nrow() > 0 and train0.ncol() > 0
    assert train1.nrow() > 0 and train1.ncol() > 0
    assert not test or (test.nrow() > 0 and test.ncol() > 0)
    assert matrixlib.are_rows_aligned(train0, train1), "matrices not aligned"
    assert not test or matrixlib.are_rows_aligned(train0, test), \
        "matrices not aligned"

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
    #return [_hash_geneid(x) for x in ids]
    # Optimization: do this without a function call.
    return [x.strip().lower() for x in ids]

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
