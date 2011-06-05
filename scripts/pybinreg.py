#!/usr/bin/env python

import os, sys

# parse_yaml_options
# update_params_from_yaml
# update_params_from_options
# 
# make_file_layout
# init_paths
# 
# read_matrices
# assert_rows_aligned
# strip_affx_control_probes
# log_matrices
# quantile_normalize_matrices
# dwd_normalize_matrices
# shiftscale_normalize_matrices
# 
# write_dataset
# write_files_for_binreg
# 
# summarize_model
# summarize_parameters
# summarize_probabilities
# summarize_signature_dataset
# summarize_signature_heatmap
# summarize_dataset_heatmap
# summarize_predictions

MAX_ANALYSES = 100   # Maximum number of analyses to do at a time.

class PyBinregParams:
    # Holds parameters for Binreg and other parameters for this
    # script.
    #
    # binreg_version       For Binreg.
    # genes                List of ints, used by Binreg one at a time.
    # metagenes            List of ints, used by Binreg one at a time.
    #
    # strip_affx           Handled by this script.
    # log_train0
    # log_train1
    # log_test
    # quantile             Handled by this script.
    # shift_scale          Handled by this script.
    # dwd                  Handled by this script.
    # dwd_bild             Handled by this script.
    # 
    # burnin               For Binreg.
    # samples              For Binreg.
    # skips                For Binreg.
    # credible_interval    For Binreg.
    # 
    # cross_validate       For Binreg.
    # noplots              Converted to make_plots for Binreg.
    # archive              Handled by this script.
    def __init__(
        self, binreg_version=2, genes=[100], metagenes=[2],
        strip_affx=False, log_train0=False, log_train1=False, log_test=False,
        quantile=False, shiftscale=False, dwd=False, dwd_bild=False,
        burnin=1000, samples=5000, skips=1, credible_interval=95,
        cross_validate=True, noplots=False, archive=False,
        ):
        self.binreg_version = binreg_version
        self.genes = genes[:]
        self.metagenes = metagenes[:]
        
        self.strip_affx = strip_affx
        self.log_train0 = log_train0
        self.log_train1 = log_train1
        self.log_test = log_test
        self.quantile = quantile
        self.shiftscale = shiftscale
        self.dwd = dwd
        self.dwd_bild = dwd_bild

        self.burnin = burnin
        self.samples = samples
        self.skips = skips
        self.credible_interval = credible_interval

        self.cross_validate = cross_validate
        self.noplots = noplots
        self.archive = archive

        self.validate()

    def validate(self):
        # Make sure the parameters are reasonable.
        # Check inputs.
        assert self.binreg_version in [1, 2]

        # Make sure genes and metagenes are in valid ranges.
        assert type(self.genes) is type([]) and len(self.genes)
        for gene in self.genes:
            assert type(gene) is type(0) and gene > 0, repr(gene)
        assert type(self.metagenes) is type([]) and len(self.metagenes)
        for metagene in self.metagenes:
            assert type(metagene) is type(0) and metagene > 0, repr(metagene)
            
        # No duplicate genes or metagenes.
        i = 0
        while i < len(self.genes):
            if self.genes[i] in self.genes[:i]:
                del self.genes[i]
            else:
                i += 1
        i = 0
        while i < len(self.metagenes):
            if self.metagenes[i] in self.metagenes[:i]:
                del self.metagenes[i]
            else:
                i += 1
        for i in range(len(self.genes)):
            assert self.genes[i] not in self.genes[i+1:]
        for i in range(len(self.metagenes)):
            assert self.metagenes[i] not in self.metagenes[i+1:]

        # Make sure the user didn't specify an unreasonable number of
        # analyses to be done.
        num_analyses = len(self.genes) * len(self.metagenes)
        if num_analyses > MAX_ANALYSES:
            assert False, "Requested %d analyses, but maximum is %d." % (
                num_analyses, MAX_ANALYSIS)

        # Make sure boolean attributes are boolean.
        bool_attrs = [
            "strip_affx", "log_train0", "log_train1", "log_test",
            "quantile", "shiftscale", "dwd", "dwd_bild",
            "cross_validate", "noplots", "archive"]
        for attr in bool_attrs:
            assert hasattr(self, attr)
            assert getattr(self, attr) in [True, False]
        
        assert self.burnin >= 0
        assert self.samples > 0
        assert self.skips >= 0
        assert self.credible_interval >= 0 and self.credible_interval <= 100

def parse_yaml_options(filename):
    import yaml
    from genomicode import filelib

    data = yaml.load(filelib.openfh(filename).read())

    # Conversion of the YAML names to the attributes for PyBinregParams.
    name2stdname = {
        "binreg_version" : "binreg_version",
        "genes" : "genes",
        "factors" : "metagenes",
        "quantile_normalize_arrays" : "quantile",
        "shift_scale" : "shiftscale",
        "burnins" : "burnin",
        "iterations" : "samples",
        "skips" : "skips",
        "remove_affymetrix_controlprobes" : "strip_affx",
        }

    # Make an object with the standard options.
    obj = filelib.GenericObject()
    for name, stdname in name2stdname.iteritems():
        assert name in data, "Missing key in YAML file: %s" % name
        assert not hasattr(obj, stdname)
        setattr(obj, stdname, data[name])

    # Make sure the values are the right type.
    obj.genes = str(obj.genes)
    obj.metagenes = str(obj.metagenes)

    znf = data.get("zero_normalized_file")
    onf = data.get("one_normalized_file")
    assert znf, "zero_normalized_file missing"
    assert onf, "one_normalized_file missing"

    # Check the data.
    assert obj.binreg_version in [1, 2]
    assert obj.quantile in [False, True]
    assert obj.shiftscale in [False, True]
    assert obj.strip_affx in [False, True]
    
    return obj, znf, onf

def update_params_from_yaml(pybr_params, filename):
    # Note: Changes pybr_params in place.
    x = parse_yaml_options(filename)
    options_yaml, train0_file, train1_file = x
    # Update the parameters based on the YAML data.
    for name in pybr_params:
        if hasattr(options_yaml, name):
            value = getattr(options_yaml, name)
            setattr(pybr_params, value)
    return pybr_params, train0_file, train1_file

def update_params_from_options(pybr_params, options, train0, train1, test):
    # test can be None.
    from genomicode import parselib
    from genomicode import binreg

    # Attributes from options to copy directly over to pybr_params.
    attributes = [
        "binreg_version",
        
        "strip_affx",
        "quantile",
        "shiftscale",
        "dwd",
        "dwd_bild",

        "burnin",
        "samples",
        "skips",
        "credible_interval",

        "noplots",
        "archive",
        ]

    for name in attributes:
        assert hasattr(pybr_params, name)
        value = getattr(options, name)
        if value is not None:
            setattr(pybr_params, name, value)

    # Parse genes and metagenes.  Accept genes as string, to allow
    # people to specify numbers like "080" for 80.
    if options.genes:
        genes = []
        for (start, end) in parselib.parse_ranges(options.genes):
            genes.extend(range(start, end+1))
        pybr_params.genes = genes
    if options.metagenes:
        metagenes = []
        for (start, end) in parselib.parse_ranges(options.metagenes):
            metagenes.extend(range(start, end+1))
        pybr_params.metagenes = metagenes

    pybr_params.log_test = False
    if options.log_the_data == "yes":
        pybr_params.log_train0 = True
        pybr_params.log_train1 = True
        if test:
            pybr_params.log_test = True
    elif options.log_the_data == "no":
        pybr_params.log_train0 = False
        pybr_params.log_train1 = False
    else:
        assert options.log_the_data == "auto"
        pybr_params.log_train0 = not binreg.is_logged_array_data(train0)
        pybr_params.log_train1 = not binreg.is_logged_array_data(train1)
        if test:
            pybr_params.log_test = not binreg.is_logged_array_data(test)

    pybr_params.validate()

    return pybr_params

def make_file_layout(outpath, num_analyses, gene, metagene):
    from genomicode import filelayout

    outpath = outpath or "."
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath = os.path.realpath(outpath)

    Path, File = filelayout.Path, filelayout.File
    
    # Have only set set of these files for the whole analysis.
    GLOBAL_FILES = [
        Path.GLOBAL_ATTIC("attic",
            File.DS_ORIG("dataset.original.gct"),
            File.DS_LOG("dataset.log.gct"),
            File.DS_QNORM("dataset.qnorm.gct"),
            File.DS_DWD("dataset.dwd.gct"),
            File.DS_DWD_BILD("dataset.dwd_bild.gct"),
            File.DS_SS("dataset.shiftscale.gct"),
            ),
        File.DS_FINAL("dataset.gct"),
        ]
    LOCAL_FILES = [
        # One set of files per set of parameters.
        Path.BINREG("binreg",
            File.BR_DESCRIPTION("description.txt"),
            File.BR_EXPRESSION("expression.txt"),
            ),
        Path.ATTIC("attic",
            File.PREDICTIONS_POV("predictions.pov"),
            File.SIGNATURE_GCT("signature.gct"),
            File.DS_SIG("dataset.sig.gct"),
            File.DS_SIG_PNG("dataset.sig.png"),
            ),
        File.PARAMETERS("parameters.txt"),
        File.PROBABILITIES("probabilities.txt"),
        File.PREDICTIONS_PNG("predictions.png"),
        File.SIGNATURE_PNG("signature.png"),
        File.MODEL("model.txt"),
        ]
    
    if num_analyses > 1:
        analysis = "%03d_GENES_%02d_MGENES" % (gene, metagene)
        
        LOCAL_FILES = [Path.ANALYSIS(analysis, *LOCAL_FILES)]
        GLOBAL_FILES = GLOBAL_FILES + [
            File.GLOBAL_PROBABILITIES("probabilities.%s.txt" % analysis),
            File.GLOBAL_PREDICTIONS_PNG("predictions.%s.png" % analysis),
            ]
    
    file_layout = Path.OUTPATH(outpath, *(GLOBAL_FILES+LOCAL_FILES))
    return file_layout

def init_paths(file_layout):
    from genomicode import filelayout

    for x in filelayout.walk(file_layout):
        dirpath, dirnames, filenames = x
        if os.path.exists(dirpath):
            continue
        os.mkdir(dirpath)

def read_matrices(train0_file, train1_file, test_file):
    from genomicode import binreg

    filenames = [train0_file, train1_file]
    if test_file:
        filenames.append(test_file)
    x = binreg.read_matrices(filenames)
    DATA, ALIGNED = x

    if not test_file:
        DATA = list(DATA) + [None]
        ALIGNED = list(ALIGNED) + [None]

    DATA_train0, DATA_train1, DATA_test = DATA
    ALIGN_train0, ALIGN_train1, ALIGN_test = ALIGNED
    
    print "train0 file has %d genes and %d samples." % DATA_train0.dim()
    print "train1 file has %d genes and %d samples." % DATA_train1.dim()
    if DATA_test:
        print "test file has %d genes and %d samples." % DATA_test.dim()
    print "Merged file has %d genes." % ALIGN_train0.nrow()
    sys.stdout.flush()

    return ALIGN_train0, ALIGN_train1, ALIGN_test

def strip_affx_control_probes(train0, train1, test):
    # test can be None
    from genomicode import binreg
    print "Stripping Affymetrix control IDs."
    train0_s = binreg.strip_affx_control_probes(train0)
    train1_s = binreg.strip_affx_control_probes(train1)
    test_s = None
    if test:
        test_s = binreg.strip_affx_control_probes(test)
    assert_rows_aligned(train0_s, train1_s, test_s)
    return train0_s, train1_s, test_s

def assert_rows_aligned(train0, train1, test):
    from genomicode import binreg
    if test:
        assert binreg.are_rows_aligned(train0, train1, test)
    else:
        assert binreg.are_rows_aligned(train0, train1)

def log_matrices(train0, train1, test, log_train0, log_train1, log_test):
    # Log each variable if necessary.  Will log in place.  Return a
    # boolean indicating whether anything was logged.  test can be None.
    from genomicode import jmath
    from genomicode import binreg

    if test is None:
        log_test = False
    
    variables = [
        ("train0", train0, log_train0),
        ("train1", train1, log_train1),
        ("test", test, log_test),
        ]
    any_files_logged = False
    for name, var, do_log in variables:
        msg = "I will not log the %s data." % name
        if do_log:
            msg = "I will log the %s data." % name
            var._X = jmath.log(var._X, base=2, safe=1)
            any_files_logged = True
        print msg
    return any_files_logged

def quantile_normalize_matrices(
    train0, train1, test, matlab=None, binreg_path=None):
    # Normalize the matrices with quantiles.  Will normalize in place.
    # test can be None.
    from genomicode import Matrix
    from genomicode import quantnorm

    assert_rows_aligned(train0, train1, test)

    X = [None] * train0.nrow()
    for i in range(len(X)):
        x = train0._X[i] + train1._X[i]
        if test:
            x = x + test._X[i]
        X[i] = x
        
    # Normalize base on the values in the training set.
    print "Normalizing with quantiles."
    I = range(train0.ncol()+train1.ncol())
    DATA = Matrix.InMemoryMatrix(X)
    # Our own quantile normalization differs from the Matlab one at
    # 4-5 decimal places.  Use the Matlab code to preserve identity.
    #DATA = quantnorm.normalize(DATA, which_columns=I)
    DATA = quantnorm.normalize_binreg(
        DATA, which_columns=I, matlab=matlab, binreg_path=binreg_path)

    # Reassign normalized values in place.
    X = DATA._X
    for i in range(len(X)):
        x = X[i]
        train0._X[i] = x[:train0.ncol()]
        train1._X[i] = x[train0.ncol():train0.ncol()+train1.ncol()]
        if test:
            test._X[i] = x[-test.ncol():]

def dwd_normalize_matrices(
    train0, train1, test, version=None, matlab=None, dwd_path=None):
    import arrayio
    from genomicode import Matrix
    from genomicode import dwdnorm

    # test must be specified to dwd normalize matrices.
    assert test, "DWD normalization requires a test set."
    
    assert_rows_aligned(train0, train1, test)
    assert train0.ncol() and train1.ncol(), "No training samples found."
    assert test.ncol(), "No test samples found."

    X = [None] * train0.nrow()
    for i in range(len(X)):
        X[i] = train0._X[i] + train1._X[i] + test._X[i]
    DATA = Matrix.InMemoryMatrix(X)

    # Y is a list of -1, 1 indicating the class of each sample.  -1 is
    # the training data, and 1 is the test data.
    Y = [-1]*(train0.ncol()+train1.ncol()) + [1]*test.ncol()
    if version == "bild":
        # Run Andrea's version.  Here, 1 is train0, -1 is train1, and
        # 2 are the tumor samples.
        Y = [1]*train0.ncol() + [-1]*train1.ncol() + [2]*test.ncol()

    if version == "bild":
        print "Normalizing with DWD (Bild)."
    else:
        print "Normalizing with DWD."
        
    DATA_n = dwdnorm.normalize(
        DATA, Y, version=version, matlab=matlab, dwd_path=dwd_path)
    assert DATA.dim() == DATA_n.dim()

    # Reassign normalized values in place.
    X_n = DATA_n._X
    for i in range(len(X_n)):
        x = X_n[i]
        train0._X[i] = x[:train0.ncol()]
        train1._X[i] = x[train0.ncol():train0.ncol()+train1.ncol()]
        test._X[i] = x[-test.ncol():]

def shiftscale_normalize_matrices(
    train0, train1, test, matlab=None, binreg_path=None):
    import arrayio
    from genomicode import Matrix
    from genomicode import shiftscalenorm

    # test must be specified to shiftscale normalize matrices.
    assert test, "Shift-scale normalization requires a test set."

    assert_rows_aligned(train0, train1, test)
    assert train0.ncol() and train1.ncol(), "No training samples found."
    assert test.ncol(), "No test samples found."

    # Shift-scale will normalize X to Y.  X is the test data and Y is
    # the training data.
    Y = [None] * train0.nrow()
    X = [None] * test.nrow()
    for i in range(len(Y)):
        Y[i] = train0._X[i] + train1._X[i]
        X[i] = test._X[i]
    DATA_X = Matrix.InMemoryMatrix(X)
    DATA_Y = Matrix.InMemoryMatrix(Y)

    print "Normalizing with shift-scale."
    DATA_n = shiftscalenorm.normalize(
        DATA_X, DATA_Y, matlab=matlab, binreg_path=binreg_path)
    assert DATA_n.dim() == DATA_X.dim()

    # Reassign normalized values in place.
    X_n = DATA_n._X
    for i in range(len(X_n)):
        test._X[i] = X_n[i]

def write_dataset(filename, train0, train1, test):
    # test can be None.
    import arrayio
    from genomicode import binreg

    matrices = [train0, train1]
    if test:
        matrices.append(test)
    DATA = binreg.merge_gct_matrices(*matrices)
    arrayio.gct_format.write(DATA, open(filename, 'w'))

def write_files_for_binreg(train0, train1, test, file_layout):
    # Format the files for binreg.  test can be None.
    from genomicode import binreg

    x = binreg.format_data_files(train0, train1, test)
    desc, exp = x
    open(file_layout.BR_DESCRIPTION, 'w').write(desc)
    open(file_layout.BR_EXPRESSION, 'w').write(exp)
    
def summarize_model(file_layout):
    filename = os.path.join(file_layout.BINREG, "genecoefficients.txt")
    assert os.path.exists(filename)

    handle = open(file_layout.MODEL, 'w')
    print >>handle, "%s\t%s" % ("Name", "Coefficient")
    for line in open(filename):
        x = line.strip().split()
        assert len(x) == 2
        coef, gene = x
        x = gene, coef
        print >>handle, "\t".join(x)
    handle.close()

def summarize_parameters(params, file_layout, num_genes, num_metagenes):
    handle = open(file_layout.PARAMETERS, 'w')
    print >>handle, "NAME\tVALUE"
    print >>handle, "Binreg Version\t%d" % params.binreg_version
    print >>handle, "Genes\t%d" % num_genes
    print >>handle, "Metagenes\t%d" % num_metagenes

    print >>handle, "Strip AFFX control\t%d" % int(params.strip_affx)
    print >>handle, "Log Train0\t%d" % int(params.log_train0)
    print >>handle, "Log Train1\t%d" % int(params.log_train1)
    print >>handle, "Log Test\t%d" % int(params.log_test)
    
    print >>handle, "Quantile Normalize\t%d" % int(params.quantile)
    print >>handle, "Shift-Scale Normalize\t%d" % int(params.shiftscale)
    print >>handle, "DWD Normalize\t%d" % int(params.dwd)
    print >>handle, "DWD Normalize (Bild)\t%d" % int(params.dwd_bild)
    
    print >>handle, "Burn In\t%d" % params.burnin
    print >>handle, "Samples\t%d" % params.samples
    print >>handle, "Skips\t%d" % params.skips
    print >>handle, "Credible Interval\t%d" % params.credible_interval

    print >>handle, "Cross Validate\t%d" % int(params.cross_validate)
    print >>handle, "Make Plots\t%d" % int(not params.noplots)
    handle.close()

def summarize_probabilities(train0, train1, test, file_layout):
    from genomicode import binreg

    x = binreg.format_predictions(
        train0, train1, test, outpath=file_layout.BINREG)
    open(file_layout.PROBABILITIES, 'w').write(x)

def summarize_signature_dataset(file_layout):
    import arrayio
    from genomicode import filelib
    from genomicode import jmath
    from genomicode import binreg

    # Read the gene_ids of interest.
    filename = os.path.join(file_layout.BINREG, "genecoefficients.txt")
    assert os.path.exists(filename)
    x = [x[1] for x in filelib.read_cols(filename)]
    gene_ids = binreg._hash_many_geneids(x)
    assert gene_ids[0] == binreg._hash_geneid("Intercept")
    gene_ids.pop(0)

    # Read the dataset.  Should be in GCT format.
    DATA = arrayio.read(file_layout.DS_FINAL)
    assert arrayio.gct_format.is_matrix(DATA)
    # Hash the IDs to make sure they match the ones in the coefficient file.
    x = binreg._hash_many_geneids(DATA.row_names("NAME"))
    DATA._row_names["NAME"] = x
    
    # Select only the signature genes.
    DATA_sig = DATA.matrix(row=gene_ids, row_header=arrayio.ROW_ID)
    assert len(gene_ids) == DATA_sig.nrow()

    # Figure out whether the data is training or test.
    type2i = { "train0" : 0, "train1" : 1, "test" : 2 }
    # Can not use sample names because they may not be unique across
    # data sets.
    #sample2type = {}
    #for d in filelib.read_row(files.probabilities, header=1):
    #    if d.Sample in sample2type:
    #        assert sample2type[d.Sample] == type2i[d.Type], \
    #               "Conflicting types for sample %s: %s %s" % (
    #            d.Sample, sample2type[d.Sample], type2i[d.Type])
    #    sample2type[d.Sample] = type2i[d.Type]
    #TYPES = [sample2type[x] for x in DATA.col_names()]
    TYPES = [None] * DATA.ncol()
    for d in filelib.read_row(file_layout.PROBABILITIES, header=1):
        d.Index = int(d.Index)
        if TYPES[d.Index] is not None:
            assert TYPES[d.Index] == type2i[d.Type], \
                   "Conflicting types for index %d: %s %s" % (
                d.Index, TYPES[d.Index], type2i[d.Type])
        TYPES[d.Index] = type2i[d.Type]

    # Pull out just the training data.
    I = [i for (i, t) in enumerate(TYPES) if t in [0, 1]]
    DATA_train = DATA_sig.matrix(None, I)
    TYPES_train = [x for x in TYPES if x in [0, 1]]

    # Sort the genes by decreasing correlation to the outcome.
    cors = [None] * DATA_train.nrow()
    for i in range(DATA_train.nrow()):
        x, y = TYPES_train, DATA_train[(i, None)]
        cors[i] = jmath.cor(x, y, safe=1)
    O = jmath.order(cors, decreasing=1)
    DATA_train = DATA_train.matrix(O, None)
    DATA_sig = DATA_sig.matrix(O, None)

    # Write out the dataset with the signature.
    handle = open(file_layout.DS_SIG, 'w')
    arrayio.gct_format.write(DATA_sig, handle)
    handle.close()

    # Normalize the data.
    #X = DATA_train._X
    #for i in range(len(X)):
    #    m = jmath.mean(X[i])
    #    # Mean center.
    #    X[i] = [x-m for x in X[i]]
    #    # Normalize to stddev of 1.
    #    s = jmath.stddev(X[i])
    #    if s == 0:
    #        continue
    #    X[i] = [x/s for x in X[i]]

    # Write it out in gct_format
    handle = open(file_layout.SIGNATURE_GCT, 'w')
    arrayio.gct_format.write(DATA_train, handle)
    handle.close()

def summarize_signature_heatmap(
    python, arrayplot, cluster, file_layout, libpath=[]):
    import arrayio
    from genomicode import plotlib

    DATA = arrayio.gct_format.read(file_layout.SIGNATURE_GCT)
    nrow, ncol = DATA.dim()

    # Make the heatmap big enough to read the gene names.  SCALE 1.2
    # is probably the minimum size, but the names are a bit fuzzy.
    SCALE = 1.5
    x = plotlib.find_tall_heatmap_size(
        nrow, ncol, min_box_height=1, min_box_width=1,
        max_box_height=40, max_box_width=40, max_total_height=768*SCALE,
        max_total_width=1024*SCALE)
    xpix, ypix = x
    #print ypix, xpix, nrow, ncol

    plotlib.plot_heatmap(
        file_layout.SIGNATURE_GCT, file_layout.SIGNATURE_PNG, xpix, ypix,
        color="bild", show_colorbar=True, show_grid=True,
        gene_center="mean", gene_normalize="var",
        gene_label=True, array_label=True,
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

def summarize_dataset_heatmap(
    python, arrayplot, cluster, file_layout, libpath=[]):
    import arrayio
    from genomicode import plotlib

    DATA = arrayio.gct_format.read(file_layout.DS_SIG)
    nrow, ncol = DATA.dim()
    
    x = plotlib.find_tall_heatmap_size(
        nrow, ncol, min_box_width=10,
        max_total_height=768*3, max_total_width=1024*3)
    xpix, ypix = x
    #print xpix, ypix, ncol, nrow

    plotlib.plot_heatmap(
        file_layout.DS_SIG, file_layout.DS_SIG_PNG, xpix, ypix,
        color="bild", show_colorbar=True, show_grid=True, 
        cluster_genes=True, gene_center="mean", gene_normalize="var",
        array_label=True, 
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

def summarize_predictions(povray, file_layout):
    # May not plot anything if there are problems with the BinReg
    # output (e.g. all predictions are nan).
    import math
    from genomicode import filelib
    from genomicode import povraygraph
    from genomicode import plotlib
    
    assert os.path.exists(file_layout.PROBABILITIES)

    COL_train0 = 0.17, 0.51, 0.99
    COL_train1 = 0.84, 0.10, 0.11
    COL_test = 0.15, 0.15, 0.16
    X, Y, error_bar = [], [], []
    pch = []
    color = []
    sample = []
    onpoint_label = []
    num_points = {}  # "test", "train0", "train1" -> count
    for d in filelib.read_row(file_layout.PROBABILITIES, header=1):
        if d.Method == "FITTED":
            continue
        sample.append(d.Sample)
        x = float(d.Metagene)
        y = float(d.Probability)

        if math.isnan(x) or math.isnan(y):
            continue
        # Calculate the error bars.
        err_l, err_u = y-float(d.Lower_CI), float(d.Upper_CI)-y
        if math.isnan(err_l):
            err_l = 0
        if math.isnan(err_u):
            err_u = 0
        # Sometimes BinReg will generate negative errors.  If it's not
        # too bad, then just ignore it.
        # Have seen negative errors as low as -0.032.
        assert err_l >= -0.05 and err_u >= -0.05, "%g %g" % (err_l, err_u)
        err_l, err_u = max(err_l, 0), max(err_u, 0)
        error = round(err_l, 3), round(err_u, 3)

        # Set the color.
        col = COL_test
        if d.Type == "train0":
            col = COL_train0
        elif d.Type == "train1":
            col = COL_train1

        # Set the shape.
        p = povraygraph.SQUARE
        if d.Type in ["train0", "train1"]:
            p = povraygraph.CIRCLE

        # Set the label.
        np = num_points.get(d.Type, 0)
        num_points[d.Type] = np+1
        onpoint_label.append(np+1)

        X.append(x)
        Y.append(y)
        pch.append(p)
        error_bar.append(error)
        color.append(col)
    assert X, "BinReg predictions didn't work.  Nothing to plot."
    assert len(X) == len(Y)
    xtick = plotlib.place_ticks(min(X), max(X))
    xtick_label = True
    ytick = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    ytick_label = ["%d%%" % int(x*100) for x in ytick]
    plot_width, plot_height = 1024, 768
    plot_width, plot_height = int(plot_width*1.5), int(plot_height*1.5)
    font = os.path.join(
        os.path.split(povraygraph.__file__)[0], "Verdana Bold.ttf")
    onpoint_label = None   # Don't display labels.
    x = povraygraph.scatter(
        X, Y, color=color, shape=pch, error_bar=error_bar, point_size=1,
        onpoint_label=onpoint_label, overpoint_label=sample,
        ylim=(0, 1.02), 
        xtick=xtick, xtick_label=xtick_label,
        ytick=ytick, ytick_label=ytick_label, tick_size=1, 
        xlabel="Metagene Score", ylabel="Probability",
        label_size=1, width=plot_width, height=plot_height,
        font=font)
    open(file_layout.PREDICTIONS_POV, 'w').write(x)
    # povray -D -J +Opredictions.png -H768 -W1024 +A0.5 predictions.pov
    r = povraygraph.povray(
        file_layout.PREDICTIONS_POV,
        outfile=file_layout.PREDICTIONS_PNG,
        height=plot_height, width=plot_width, antialias=0.5, quality=9,
        povray_bin=povray)
    print r.read()
    assert os.path.exists(file_layout.PREDICTIONS_PNG), \
           "Failed to plot predictions."

## def make_FILES(files):
##     handle = open(files.FILES, 'w')
##     w = handle.write
##     w("This describes the files created by the analysis.\n")
##     w("\n")
##     w("Summary of results:\n")
##     w("probabilities.txt     Predicted probabilities on all samples.\n")
##     w("predictions.png       Plot of predictions.\n")
##     w("signature.png         Heatmap of the signature.\n")
##     w("FILES                 This file.\n")
    
##     w("\n")
##     w("Raw results:\n")
##     w("normalizedX.txt       Normalized expression data.\n")
##     w("normalizedX.pcl       Normalized expression data in PCL format.\n")
##     w("signature.pcl         Expression profile of the signature in PCL "
##       "format.\n")
##     w("signature.cdt         Same as signature.pcl.\n")
##     w("trainingcases.txt     Predictions on training set (fitted).\n")
##     w("crossvalidation.txt   Predictions on training set "
##       "(cross-validation).\n")
##     w("validationcases.txt   Predictions on test set.\n")
##     w("genecoefficients.txt  Parameters for regression model.\n")
##     w("plotdata.mat          Internal data for generating plots.\n")
##     w("predictions.pov       POV-Ray description file.\n")
##     w("\n")
##     w("Input to BinReg:\n")
##     w("description.txt       Gene IDs.\n")
##     w("expression.txt        Expression data.\n")
##     w("preferences.dat       BinReg preference file.\n")
##     handle.close()

def main():
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options] train0_file train1_file [test_file]"
    parser = OptionParser(usage=usage, version="%prog 01")

    pybr_params = PyBinregParams()

    parser.add_option(
        "-v", "--binreg_version", dest="binreg_version", type="int",
        default=None,
        help="Use binreg version 1 or 2 (default %d)." %
        pybr_params.binreg_version)
    parser.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python (optional).")
    parser.add_option(
        "", "--dwd_path", dest="dwd_path", default=None,
        help="Specify the path of BatchAdjust.")
    parser.add_option(
        "", "--binreg", dest="binreg_path", default=None,
        help="Specify the path to the BinReg2.0 code.")
    parser.add_option(
        "", "--matlab", dest="matlab", default="matlab",
        help="Specify the command to run matlab.")
    parser.add_option(
        "", "--arrayplot", dest="arrayplot", default="arrayplot.py",
        help="Specify the command to run arrayplot.")
    parser.add_option(
        "", "--povray", dest="povray", default="povray",
        help="Specify the command to run povray.")
    parser.add_option(
        "", "--cluster", dest="cluster", default=None,
        help="Specify the command to run cluster.")
    parser.add_option(
        "", "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_option(
        "-o", "--outpath", dest="outpath", type="string", default=None,
        help="Save files in this path.")
    parser.add_option(
        "-z", "", dest="archive", action="store_true", default=None,
        help="Archive the raw output.  Helpful for GenePattern.")
    parser.add_option(
        "", "--noplots", dest="noplots", action="store_true", default=None,
        help="Do not make any plots.")
    parser.add_option(
        "", "--yaml", dest="yaml", default=None,
        help="Configure the analysis from a YAML-formatted file.  "
        "Command line parameters takes precedence over the ones from "
        "the file.")
    
    group = OptionGroup(parser, "Signature")
    group.add_option(
        "-g", "--genes", dest="genes", default=None,
        help="Number of genes to use (default %s)." % pybr_params.genes[0])
    group.add_option(
        "-m", "--metagenes", dest="metagenes", default=None,
        help="Number of metagenes to use (default %s)." %
        pybr_params.metagenes[0])
    parser.add_option_group(group)

    group = OptionGroup(parser, "Normalization")
    group.add_option(
        "", "--strip_affx", dest="strip_affx", action="store_true",
        default=None,
        help="Strip the AFFX control probes.")
    group.add_option(
        "-l", "--log_the_data", dest="log_the_data", type="choice",
        choices=["yes", "no", "auto"], default="auto",
        help="Log the data before analyzing.  "
        "Must be 'yes', 'no', or 'auto' (default).")
    group.add_option(
        "-q", "--quantile", dest="quantile", action="store_true",
        default=None,
        help="Quantile normalize the data.")
    group.add_option(
        "-d", "--dwd", dest="dwd", action="store_true",
        default=None,
        help="DWD normalize the data.")
    group.add_option(
        "", "--dwd_bild", dest="dwd_bild", action="store_true",
        default=None,
        help="DWD normalize the data, using the Bild method.")
    group.add_option(
        "-s", "--shiftscale", dest="shiftscale", action="store_true",
        default=None,
        help="Shift-scale normalize the data.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "MCMC Internals")
    group.add_option(
        "", "--burnin", dest="burnin", type="int", default=None,
        help="Number of burn-in samples (default %d)." % pybr_params.burnin)
    group.add_option(
        "", "--samples", dest="samples", type="int", default=None,
        help="Number of samples (default %d)." % pybr_params.samples)
    group.add_option(
        "", "--skips", dest="skips", type="int", default=None,
        help="Number of skips (default %d)." % pybr_params.skips)
    group.add_option(
        "", "--credible_interval", dest="credible_interval", type="int",
        default=None,
        help="Credible interval (default %d)." % pybr_params.credible_interval)
    parser.add_option_group(group)

    # Parse the arguments.
    options, args = parser.parse_args()

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import after the library path is set.
    import arrayio
    from genomicode import filelib
    from genomicode import archive
    from genomicode import binreg
    from genomicode import genepattern

    genepattern.fix_environ_path()

    # If YAML file is specified, then use the values stored in the
    # YAML file as default.  Precedence is USER > YAML > DEFAULTS.
    if options.yaml:
        # Since train0_file and train1_file are specified in the YAML
        # file, args should only contain test_file.
        if len(args) != 1:
            parser.error("Please specify test file.")
        test_file, = args
        x = update_params_from_yaml(pybr_params, options.yaml)
        pybr_params, train0_file, train1_file = x
        # Set args to the training and test files.
        args = train0_file, train1_file, test_file

    if len(args) == 2:
        train0_file, train1_file = args
        test_file = None
    elif len(args) == 3:
        train0_file, train1_file, test_file = args
    elif len(args) < 2:
        parser.error("Please specify train0 and train1 files.")
    else:
        parser.error("Too many files.")
        
    assert filelib.exists(train0_file), "File not found: %s" % train0_file
    assert filelib.exists(train1_file), "File not found: %s" % train1_file
    if test_file:
        assert filelib.exists(test_file), "File not found: %s" % test_file

    if not test_file:
        if options.shiftscale:
            print "No test file, so disabling shift-scale normalization."
            options.shiftscale = False
        if options.dwd:
            print "No test file, so disabling DWD normalization."
            options.dwd = False
        if options.dwd_bild:
            print "No test file, so disabling DWD (Bild) normalization."
            options.dwd_bild = False

    # Read each of the input files and align them.
    x = read_matrices(train0_file, train1_file, test_file)
    train0, train1, test = x

    x = update_params_from_options(pybr_params, options, train0, train1, test)
    pybr_params = x

    num_analyses = len(pybr_params.genes) * len(pybr_params.metagenes)
    file_layout = make_file_layout(
        options.outpath, num_analyses,
        pybr_params.genes[0], pybr_params.metagenes[0])
    init_paths(file_layout)

    # To debug creation of scatterplot.
    #summarize_predictions(options.povray, file_layout); return

    if pybr_params.strip_affx:
        x = strip_affx_control_probes(train0, train1, test)
        train0, train1, test = x

    # Standardize each of the matrices to GCT format.
    train0 = arrayio.convert(train0, to_format=arrayio.gct_format)
    train1 = arrayio.convert(train1, to_format=arrayio.gct_format)
    if test:
        test = arrayio.convert(test, to_format=arrayio.gct_format)
    write_dataset(file_layout.DS_ORIG, train0, train1, test)

    # Log normalize each of the files, if necessary.  Will log in place.
    logged = log_matrices(
        train0, train1, test, pybr_params.log_train0, pybr_params.log_train1,
        pybr_params.log_test)
    if logged:
        write_dataset(file_layout.DS_LOG, train0, train1, test)

    if pybr_params.quantile:
        quantile_normalize_matrices(
            train0, train1, test, matlab=options.matlab,
            binreg_path=options.binreg_path)
        write_dataset(file_layout.DS_QNORM, train0, train1, test)

    if pybr_params.dwd:
        dwd_normalize_matrices(
            train0, train1, test, matlab=options.matlab,
            dwd_path=options.dwd_path)
        write_dataset(file_layout.DS_DWD, train0, train1, test)

    if pybr_params.dwd_bild:
        dwd_normalize_matrices(
            train0, train1, test, version="bild", matlab=options.matlab,
            dwd_path=options.dwd_path)
        write_dataset(file_layout.DS_DWD_BILD, train0, train1, test)

    if pybr_params.shiftscale:
        shiftscale_normalize_matrices(
            train0, train1, test, matlab=options.matlab,
            binreg_path=options.binreg_path)
        write_dataset(file_layout.DS_SS, train0, train1, test)

    print "Writing final dataset."
    write_dataset(file_layout.DS_FINAL, train0, train1, test)

    # Make a list of the analyses to do.
    analyses = []  # list of (gene, metagene)
    for gene in pybr_params.genes:
        assert gene <= train0.nrow(), \
            "Merged data set has %d genes, but model requires %d." % (
            train0.nrow(), gene)
        for metagene in pybr_params.metagenes:
            assert metagene <= gene
            analyses.append((gene, metagene))

    # Now do each analysis.
    for gene, metagene in analyses:
        x = "Analyzing with %d genes and %d metagenes." % (gene, metagene)
        if metagene == 1:
            x = x.replace("metagenes", "metagene")
        print x
        
        # Generate the file layout for this analysis.
        file_layout = make_file_layout(
            options.outpath, num_analyses, gene, metagene)
        init_paths(file_layout)

        # Format the parameters and output files for binreg.
        write_files_for_binreg(train0, train1, test, file_layout)
        params = binreg.BinregParams(
            binreg_version=pybr_params.binreg_version,
            cross_validate=int(pybr_params.cross_validate),
            make_plots=0,
            num_genes=gene, num_metagenes=metagene,
            quantile_normalize=0, shift_scale_normalize=0,
            num_burnin=pybr_params.burnin, num_iterations=pybr_params.samples,
            num_skips=pybr_params.skips,
            credible_interval=pybr_params.credible_interval)

        # Run Binreg.
        r = binreg.binreg_raw(
            file_layout.BR_EXPRESSION, file_layout.BR_DESCRIPTION,
            1, params, matlab=options.matlab,
            binreg_path=options.binreg_path, outpath=file_layout.BINREG)
        print r.read()
        sys.stdout.flush()

        # Generate some files for output.
        summarize_model(file_layout)
        summarize_parameters(pybr_params, file_layout, gene, metagene)
        summarize_probabilities(train0, train1, test, file_layout)
        summarize_signature_dataset(file_layout)
        
        # Make some plots, if desired.
        if not pybr_params.noplots:
            summarize_signature_heatmap(
                options.python, options.arrayplot, options.cluster,
                file_layout, options.libpath)
            summarize_dataset_heatmap(
                options.python, options.arrayplot, options.cluster,
                file_layout, options.libpath)
            summarize_predictions(options.povray, file_layout)
        sys.stdout.flush()

        # Don't need, because archive takes care of it.
        ## Compress the files if they are too big.
        #pybr_params.compression_limit = 10  # gzip if over 10 Mb.
        #if pybr_params.compression_limit:
        #    compress(pybr_params.compression_limit, file_layout)

        if num_analyses <= 1:
            continue
        # Now do some cleanup if multiple analyses were requested.

        # If there were multiple genes and metagenes specified,
        # make a copy of some files for convenience.  Need to copy
        # it rather than symlink, in case we archive it later.
        fl = file_layout
        files_to_copy = [
            (fl.PROBABILITIES, fl.GLOBAL_PROBABILITIES),
            (fl.PREDICTIONS_PNG, fl.GLOBAL_PREDICTIONS_PNG),
            ]
        for src, dst in files_to_copy:
            assert os.path.exists(src)
            os.system("cp -p '%s' '%s'" % (src, dst))

        if pybr_params.archive:
            archive.zip_path(file_layout.ANALYSIS)
        sys.stdout.flush()

    if pybr_params.archive:
        if os.path.exists(file_layout.BINREG):
            archive.zip_path(file_layout.BINREG)
        if os.path.exists(file_layout.ATTIC):
            archive.zip_path(file_layout.ATTIC)
        if os.path.exists(file_layout.GLOBAL_ATTIC):
            archive.zip_path(file_layout.GLOBAL_ATTIC)

    #make_FILES(files)
    print "Done."
    
if __name__ == '__main__':
    main()
    #import profile; profile.run("main()")
