#!/usr/bin/env python

import os, sys

# init_params
# make_params_from_default
# make_params_from_yaml
# make_params_from_options
# check_params
# 
# make_analyses
#
# make_file_layout
# init_paths
# 
# read_matrices
# read_matrices_all
# assert_rows_aligned
# strip_affx_control_probes
# log_matrices
# quantile_normalize_matrices
# dwd_normalize_matrices
# shiftscale_normalize_matrices
# _shiftscale_normalize_standard
# _shiftscale_normalize_normref
# 
# write_dataset
# write_files_for_binreg
# run_analysis
# _run_binreg_standard
# _run_binreg_normref
# 
# summarize_model
# summarize_parameters
# summarize_probabilities
# summarize_signature_dataset
# summarize_signature_heatmap
# summarize_dataset_heatmap
# summarize_predictions
# summarize_report
# summarize_all_report
# _make_signature_figure
# _make_prediction_figure
# _make_footer
#
# pretty_runtime
# pretty_hostname

MAX_ANALYSES = 128   # Maximum number of analyses to do at a time.

class Analysis:
    def __init__(self, params, subpath):
        self.params = params
        # subpath should be None if this analysis should be in the
        # root directory.
        self.subpath = subpath

class PyBinregParams:
    # Holds parameters for Binreg and other parameters for this
    # script.
    #
    # binreg_version       For Binreg.
    # genes                List of ints, used by Binreg one at a time.
    # metagenes            List of ints, used by Binreg one at a time.
    #
    # train0_file
    # train1_file
    # test_file
    # 
    # log_train0
    # log_train1
    # log_test
    # log_normref
    # strip_affx           Handled by this script.
    # quantile             Handled by this script.
    # shift_scale          Handled by this script.
    # dwd                  Handled by this script.
    # dwd_bild             Handled by this script.
    # normref_file         Handled by this script.
    #
    # burnin               For Binreg.
    # samples              For Binreg.
    # skips                For Binreg.
    # credible_interval    For Binreg.
    # 
    # cross_validate       For Binreg.
    # noplots              Converted to make_plots for Binreg.
    # label_samples        Handled by this script.
    # draw_error_bars      Handled by this script.
    # archive              Handled by this script.
    def __init__(
        self, binreg_version=None, genes=None, metagenes=None,
        train0_file=None, train1_file=None, test_file=None,
        normref_file=None, strip_affx=None, 
        log_train0=None, log_train1=None, log_test=None, log_normref=None,
        quantile=None, shiftscale=None, dwd=None, dwd_bild=None,
        burnin=None, samples=None, skips=None, credible_interval=None,
        cross_validate=None, noplots=None, label_samples=None,
        draw_error_bars=None, archive=None,
        ):
        self.binreg_version = binreg_version
        self.genes = genes
        self.metagenes = metagenes

        self.train0_file = train0_file
        self.train1_file = train1_file
        self.test_file = test_file
        self.normref_file = normref_file
        
        self.strip_affx = strip_affx
        self.log_train0 = log_train0
        self.log_train1 = log_train1
        self.log_test = log_test
        self.log_normref = log_normref
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
        self.label_samples = label_samples
        self.draw_error_bars = draw_error_bars
        self.archive = archive

def init_params():
    # Initialize the parameters.
    from optparse import OptionParser, OptionGroup
    
    # Initialize the default parameters.  Needed for some of the
    # options.
    default_params = make_params_from_default()

    usage = "usage: %prog [options] train0_file train1_file [test_file]"
    parser = OptionParser(usage=usage, version="%prog 02")

    parser.add_option(
        "-v", "--binreg_version", dest="binreg_version", type="int",
        default=None,
        help="Use binreg version 1 or 2 (default %d)." %
        default_params.binreg_version)
    parser.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python (optional).")
    parser.add_option(
        "", "--matlab", dest="matlab", default="matlab",
        help="Specify the command to run matlab.")
    parser.add_option(
        "", "--povray", dest="povray", default="povray",
        help="Specify the command to run povray.")
    parser.add_option(
        "", "--cluster", dest="cluster", default=None,
        help="Specify the command to run cluster.")
    parser.add_option(
        "", "--dwd_path", dest="dwd_path", default=None,
        help="Specify the path of BatchAdjust.")
    parser.add_option(
        "", "--binreg", dest="binreg_path", default=None,
        help="Specify the path to the BinReg2.0 code.")
    parser.add_option(
        "", "--arrayplot", dest="arrayplot", default=None,
        help="Specify the command to run arrayplot.")
    parser.add_option(
        "", "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_option(
        "-o", "--outpath", dest="outpath", type="string", default=None,
        help="Save files in this path.")
    parser.add_option(
        "-z", "", dest="archive", action="store_true", default=None,
        help="Archive the raw output.  Helpful for GenePattern.")
    #parser.add_option(
    #    "-j", "", dest="num_procs", type="int", default=1,
    #    help="Number of jobs to run in parallel.")
    parser.add_option(
        "", "--yaml", dest="yaml", default=None,
        help="Configure the analysis from a YAML-formatted file.  "
        "Command line parameters takes precedence over the ones from "
        "the file.")
    
    group = OptionGroup(parser, "Plotting")
    group.add_option(
        "", "--noplots", dest="noplots", action="store_true", default=None,
        help="Do not make any plots.")
    group.add_option(
        "", "--label_samples", dest="label_samples", type="choice",
        choices=["yes", "no", "auto"], default="auto",
        help="Label the samples in the scatter plot.  "
        "Must be 'yes', 'no', or 'auto' (default).")
    group.add_option(
        "", "--draw_error_bars", dest="draw_error_bars", type="choice",
        choices=["yes", "no", "auto"], default="auto",
        help="Draw the error bars in the scatter plot.  "
        "Must be 'yes', 'no', or 'auto' (default).")
    group.add_option(
        "", "--credible_interval", dest="credible_interval", type="int",
        default=None, help="Credible interval (default %d)." %
        default_params.credible_interval)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Signature")
    group.add_option(
        "-g", "--genes", dest="genes", default=None,
        help="Number of genes to use (default %s)." % default_params.genes)
    group.add_option(
        "-m", "--metagenes", dest="metagenes", default=None,
        help="Number of metagenes to use (default %s)." %
        default_params.metagenes)
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
    group.add_option(
        "--normalization_reference_file", dest="normref_file",
        type="string", default=None,
        help="Use this file to help normalize the test file.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "MCMC Internals")
    group.add_option(
        "", "--burnin", dest="burnin", type="int", default=None,
        help="Number of burn-in samples (default %d)." % default_params.burnin)
    group.add_option(
        "", "--samples", dest="samples", type="int", default=None,
        help="Number of samples (default %d)." % default_params.samples)
    group.add_option(
        "", "--skips", dest="skips", type="int", default=None,
        help="Number of skips (default %d)." % default_params.skips)
    parser.add_option_group(group)

    # Parse the arguments.
    options, args = parser.parse_args()

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import after the library path is set.
    from genomicode import genepattern

    genepattern.fix_environ_path()

    yaml_params = None
    if options.yaml:
        yaml_params = make_params_from_yaml(options.yaml, args)
        # HACK: If YAML, then args will only consist of the test_file.
        # The train0_file and train1_file will be provided in the YAML
        # file.  However, the options parameters might depend on
        # these.  If this is the case, then add the train0_file and
        # train1_file to args.
        assert yaml_params.train0_file and yaml_params.train_file and \
               yaml_params.test_file
        x = yaml_params.train0_file, yaml_params.train1_file, \
            yaml_params.test_file
        args = x
    option_params_list = make_params_from_options(options, args)

    # Precedence is OPTION > YAML > DEFAULT.
    param_list = []
    for option_params in option_params_list:
        x = [option_params, yaml_params, default_params]
        parameters = [x for x in x if x]
        keywds = {}
        for i in range(len(parameters)-1, -1, -1):
            p = parameters[i]
            for name in dir(p):
                if name.startswith("_"):
                    continue
                value = getattr(p, name)
                if value is None:
                    continue
                keywds[name] = value
        x = PyBinregParams(**keywds)
        # Make sure params are complete and valid.
        check_params(x)
        param_list.append(x)

    # Make sure the user didn't specify an unreasonable number of
    # analyses to be done.
    if len(param_list) > MAX_ANALYSES:
        assert False, "Requested %d analyses, but maximum is %d." % (
            len(param_list), MAX_ANALYSIS)
    
    return options, param_list

def make_params_from_default():
    x = PyBinregParams(
        binreg_version=2,
        genes=100,
        metagenes=2,
        normref_file=None,
        strip_affx=False,
        log_train0=False,
        log_train1=False,
        log_test=False,
        log_normref=False,
        quantile=False,
        shiftscale=False,
        dwd=False,
        dwd_bild=False,
        burnin=1000,
        samples= 5000,
        skips=1,
        credible_interval=95,
        cross_validate=True,
        noplots=False,
        label_samples=True,
        draw_error_bars=True,
        archive=False,
        )
    return x

def make_params_from_yaml(filename, args):
    import yaml
    from genomicode import filelib

    # Since train0_file and train1_file are specified in the YAML
    # file, args should only contain test_file.
    assert len(args) == 1, "Please specify test file."
    test_file, = args

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
    keywds = {}
    for name, stdname in name2stdname.iteritems():
        assert name in data, "Missing key in YAML file: %s" % name
        assert stdname not in keywds
        keywds[stdname] = data[name]

    # Make sure the values are the right type.
    keywds["genes"] = int(obj.genes)
    keywds["metagenes"] = int(obj.metagenes)

    znf = data.get("zero_normalized_file")
    onf = data.get("one_normalized_file")
    assert znf, "zero_normalized_file missing"
    assert onf, "one_normalized_file missing"
    keywds["train0_file"] = znf
    keywds["train1_file"] = onf
    keywds["test_file"] = test_file

    # Check the data.
    assert keywds["binreg_version"] in [1, 2]
    assert keywds["quantile"] in [False, True]
    assert keywds["shiftscale"] in [False, True]
    assert keywds["strip_affx"] in [False, True]

    x = PyBinregParams(**keywds)
    return x

def make_params_from_options(options, args):
    # Return a list of PyBinregParams.
    import copy
    import itertools
    from genomicode import parselib
    from genomicode import binreg

    options = copy.copy(options)
    
    # Parse the file arguments.
    train0_file = train1_file = test_file = None
    if len(args) == 2:
        train0_file, train1_file = args
    elif len(args) == 3:
        train0_file, train1_file, test_file = args
    elif len(args) < 2:
        raise AssertionError, "Please specify train0 and train1 files."
    else:
        raise AssertionError, "Too many files."
    normref_file = options.normref_file

    # Parse genes and metagenes.  Accept genes as string, to allow
    # people to specify numbers like "080" for 80.
    if options.genes:
        genes = []
        for (start, end) in parselib.parse_ranges(options.genes):
            genes.extend(range(start, end+1))
    if options.metagenes:
        metagenes = []
        for (start, end) in parselib.parse_ranges(options.metagenes):
            metagenes.extend(range(start, end+1))

    # Normalization reference file only works with shift-scale.
    if options.normref_file and not options.shiftscale:
        if test_file:
            print("Normalization reference provided, so enabling shift-scale "
                  "normalization.")
            options.shiftscale = True

    # Configure the normalization options based on the test file.
    if not test_file:
        if options.shiftscale:
            print "No test file, so disabling shift-scale normalization."
            options.shiftscale = False
            options.normref_file = None
        if options.dwd:
            print "No test file, so disabling DWD normalization."
            options.dwd = False
        if options.dwd_bild:
            print "No test file, so disabling DWD (Bild) normalization."
            options.dwd_bild = False
            
    # Read each of the input files and align them.
    x = read_matrices(train0_file, train1_file, test_file, normref_file)
    train0, train1, test, normref = x

    # Configure the logging options based on the files.
    log_test = log_normref = False
    if options.log_the_data == "yes":
        log_train0 = True
        log_train1 = True
        if test:
            log_test = True
        if normref:
            log_normref = True
    elif options.log_the_data == "no":
        log_train0 = False
        log_train1 = False
    else:
        assert options.log_the_data == "auto"
        log_train0 = not binreg.is_logged_array_data(train0)
        log_train1 = not binreg.is_logged_array_data(train1)
        if test:
            log_test = not binreg.is_logged_array_data(test)
        if normref:
            log_normref = not binreg.is_logged_array_data(normref)

    # Configure the plotting parameters.
    assert options.label_samples in ["yes", "no", "auto"]
    assert options.draw_error_bars in ["yes", "no", "auto"]
    label_samples = True
    if options.label_samples == "no" or (
        options.label_samples == "auto" and test and test.ncol() > 50):
        # Don't label if it's too crowded.
        label_samples = False
    draw_error_bars = True
    if options.draw_error_bars == "no":
        draw_error_bars = False

    # Set some attributes that can be directly copied from the options.
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
    keywds = {}
    keywds["train0_file"] = train0_file
    keywds["train1_file"] = train1_file
    keywds["test_file"] = test_file
    keywds["normref_file"] = normref_file
    keywds["log_train0"] = log_train0
    keywds["log_train1"] = log_train1
    keywds["log_test"] = log_test
    keywds["log_normref"] = log_normref
    keywds["label_samples"] = label_samples
    keywds["draw_error_bars"] = draw_error_bars
    for name in attributes:
        assert name not in keywds
        value = getattr(options, name)
        if value is not None:
            keywds[name] = value

    paramlist = []
    seen = {}   # No duplicate genes or metagenes.
    for x in itertools.product(genes, metagenes):
        g, mg = x
        if (g, mg) in seen:
            continue
        seen[(g, mg)] = 1
        
        keywds["genes"] = g
        keywds["metagenes"] = mg
        x = PyBinregParams(**keywds)
        paramlist.append(x)
    return paramlist

def check_params(params):
    # Make sure the parameters are reasonable.
    from genomicode import filelib
    
    assert params.binreg_version in [1, 2], repr(params.binreg_version)

    # Make sure genes and metagenes are in valid ranges.
    assert type(params.genes) is type(0) and params.genes > 0, \
           repr(params.genes)
    assert type(params.metagenes) is type(0) and params.metagenes > 0, \
           repr(params.metagenes)

    # Make sure each file exists.
    assert filelib.exists_nz(params.train0_file), \
           "File not found: %s" % params.train0_file
    assert filelib.exists_nz(params.train1_file), \
           "File not found: %s" % params.train1_file
    if params.test_file:
        assert filelib.exists_nz(params.test_file), \
               "File not found: %s" % params.test_file
    if params.normref_file:
        assert filelib.exists_nz(params.normref_file), \
               "File not found: %s" % params.normref_file
        
    # Make sure boolean attributes are boolean.
    bool_attrs = [
        "strip_affx", "log_train0", "log_train1", "log_test", "log_normref",
        "quantile", "shiftscale", "dwd", "dwd_bild",
        "cross_validate", "noplots", "label_samples", "draw_error_bars",
        "archive"]
    for attr in bool_attrs:
        assert hasattr(params, attr)
        assert getattr(params, attr) in [True, False]

    assert params.burnin >= 0
    assert params.samples > 0
    assert params.skips >= 0
    assert params.credible_interval >= 0 and params.credible_interval <= 100

def make_analyses(param_list):
    # Return a list of Analysis objects.
    analyses = []
    for params in param_list:
        path = None
        if len(param_list) > 1:
            path = "%03d_GENES_%02d_MGENES" % (
                params.genes, params.metagenes)
        x = Analysis(params, path)
        analyses.append(x)
            
    return analyses

def make_file_layout(outpath, analysis_path):
    from genomicode import filelayout

    outpath = outpath or "."
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath = os.path.realpath(outpath)

    Path, File = filelayout.Path, filelayout.File
    
    # Have only one set of these files for the whole analysis.
    GLOBAL_FILES = [
        Path.GLOBAL_ATTIC("attic",
            File.DS_ORIG("dataset.original.gct"),
            File.DS_LOG("dataset.log.gct"),
            File.DS_QNORM("dataset.qnorm.gct"),
            File.DS_DWD("dataset.dwd.gct"),
            File.DS_DWD_BILD("dataset.dwd_bild.gct"),
            File.DS_SS("dataset.shiftscale.gct"),
            File.NR_ORIG("normref.original.gct"),
            File.NR_QNORM("normref.qnorm.gct"),
            File.NR_FINAL("normref.final.gct"),
            ),
        File.DS_FINAL("dataset.gct"),
        ]
    LOCAL_FILES = [
        # One set of files per set of parameters.
        Path.BINREG("binreg",
            File.BR_DESCRIPTION("description.txt"),
            File.BR_EXPRESSION("expression.txt"),
            File.BR_COEFFICIENTS("genecoefficients.txt"),
            File.BR_VALIDATION("validationcases.txt"),
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
        File.REPORT("REPORT.html"),
        ]

    if analysis_path:
        LOCAL_FILES = [Path.ANALYSIS(analysis_path, *LOCAL_FILES)]
        GLOBAL_FILES = GLOBAL_FILES + [
            File.GLOBAL_REPORT("REPORT.html"),
            File.GLOBAL_PROBABILITIES("probabilities.%s.txt" % analysis_path),
            File.GLOBAL_PREDICTIONS_PNG("predictions.%s.png" % analysis_path),
            File.GLOBAL_SIGNATURE_PNG("signature.%s.png" % analysis_path),
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

MATRIX_CACHE = {}
def read_matrices_all(train0_file, train1_file, test_file, normref_file):
    global MATRIX_CACHE
    from genomicode import parselib
    from genomicode import binreg

    filenames = [train0_file, train1_file]
    if test_file:
        filenames.append(test_file)
    if normref_file:
        filenames.append(normref_file)
    x = binreg.read_matrices(filenames, cache=MATRIX_CACHE)
    DATA, ALIGNED = x

    DATA_train0, DATA_train1 = DATA[:2]
    ALIGN_train0, ALIGN_train1 = ALIGNED[:2]
    DATA_test = ALIGN_test = DATA_normref = ALIGN_normref = None
    if test_file:
        DATA_test = DATA[2]
        ALIGN_test = ALIGNED[2]
    if normref_file:
        DATA_normref = DATA[-1]
        ALIGN_normref = ALIGNED[-1]

    x = (DATA_train0, DATA_train1, DATA_test, DATA_normref,
         ALIGN_train0, ALIGN_train1, ALIGN_test, ALIGN_normref)
    return x

def read_matrices(train0_file, train1_file, test_file, normref_file):
    x = read_matrices_all(train0_file, train1_file, test_file, normref_file)
    (DATA_train0, DATA_train1, DATA_test, DATA_normref,
     ALIGN_train0, ALIGN_train1, ALIGN_test, ALIGN_normref) = x
    return ALIGN_train0, ALIGN_train1, ALIGN_test, ALIGN_normref

def assert_rows_aligned(train0, train1, test, normref):
    from genomicode import binreg
    x = [train0, train1, test, normref]
    matrices = [x for x in x if x]
    assert binreg.are_rows_aligned(*matrices)

def strip_affx_control_probes(train0, train1, test, normref):
    # test, normref can be None
    from genomicode import binreg
    
    print "Stripping Affymetrix control IDs."
    train0_s = binreg.strip_affx_control_probes(train0)
    train1_s = binreg.strip_affx_control_probes(train1)
    test_s = None
    if test:
        test_s = binreg.strip_affx_control_probes(test)
    normref_s = None
    if normref:
        normref_s = binreg.strip_affx_control_probes(normref)
    assert_rows_aligned(train0_s, train1_s, test_s, normref_s)
    return train0_s, train1_s, test_s, normref_s

def log_matrices(
    train0, train1, test, normref,
    log_train0, log_train1, log_test, log_normref):
    # Log each variable if necessary.  Will log in place.  Return a
    # boolean indicating whether anything was logged.  test can be None.
    from genomicode import jmath
    from genomicode import binreg

    variables = [
        ("train0", train0, log_train0),
        ("train1", train1, log_train1),
        ("test", test, log_test),
        ("normref", normref, log_normref),
        ]
    any_files_logged = False
    for name, var, do_log in variables:
        if var is None:
            continue
        msg = "I will not log the %s data." % name
        if do_log:
            msg = "I will log the %s data." % name
            var._X = jmath.log(var._X, base=2, safe=1)
            any_files_logged = True
        print msg
    sys.stdout.flush()
    return any_files_logged

def quantile_normalize_matrices(
    train0, train1, test, normref, matlab=None, binreg_path=None):
    # Normalize the matrices with quantiles.  Will normalize in place.
    # test can be None.
    from genomicode import Matrix
    from genomicode import quantnorm

    assert_rows_aligned(train0, train1, test, None)

    X = [None] * train0.nrow()
    for i in range(len(X)):
        x = train0._X[i] + train1._X[i]
        if test:
            x = x + test._X[i]
        if normref:
            x = x + normref._X[i]
        X[i] = x
        
    # Normalize base on the values in the training set.
    print "Normalizing with quantiles."
    I = range(train0.ncol()+train1.ncol())
    DATA = Matrix.InMemoryMatrix(X)
    # Our own quantile normalization differs from the Matlab one at
    # the 4-5th decimal place.  Use the Matlab code so that the
    # results will be exactly the same.
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
            j = train0.ncol()+train1.ncol()
            test._X[i] = x[j:j+test.ncol()]
        if normref:
            j = train0.ncol()+train1.ncol()
            if test:
                j += test.ncol()
            normref._X[i] = x[j:j+normref.ncol()]

def dwd_normalize_matrices(
    train0, train1, test, version=None, matlab=None, dwd_path=None):
    import arrayio
    from genomicode import Matrix
    from genomicode import dwdnorm

    # test must be specified to dwd normalize matrices.
    assert test, "DWD normalization requires a test set."
    
    assert_rows_aligned(train0, train1, test, None)
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
    train0, train1, test, normref, matlab=None, binreg_path=None):
    if not normref:
        _shiftscale_normalize_standard(
            train0, train1, test, matlab=matlab, binreg_path=binreg_path)
    else:
        _shiftscale_normalize_normref(
            train0, train1, test, normref, matlab=matlab,
            binreg_path=binreg_path)

def _shiftscale_normalize_standard(
    train0, train1, test, matlab=None, binreg_path=None):
    import arrayio
    from genomicode import Matrix
    from genomicode import shiftscalenorm

    print "Normalizing with shift-scale."; sys.stdout.flush()

    # test must be specified to shiftscale normalize matrices.
    assert test, "Shift-scale normalization requires a test set."
    assert train0.ncol() and train1.ncol(), "Training samples missing."
    assert test.ncol(), "No test samples found."
    assert_rows_aligned(train0, train1, test, None)

    # Shift-scale will normalize X to Y.  X is the test data and Y is
    # the training data.
    Y = [None] * train0.nrow()
    X = [None] * test.nrow()
    for i in range(len(Y)):
        Y[i] = train0._X[i] + train1._X[i]
        X[i] = test._X[i]
    DATA_X = Matrix.InMemoryMatrix(X)
    DATA_Y = Matrix.InMemoryMatrix(Y)

    DATA_n = shiftscalenorm.normalize(
        DATA_X, DATA_Y, matlab=matlab, binreg_path=binreg_path)
    assert DATA_n.dim() == DATA_X.dim()

    # Reassign normalized values in place.
    X_n = DATA_n._X
    for i in range(len(X_n)):
        test._X[i] = X_n[i]

def _shiftscale_normalize_normref(
    train0, train1, test, normref, matlab=None, binreg_path=None):
    import arrayio
    from genomicode import Matrix
    from genomicode import shiftscalenorm

    # test must be specified to shiftscale normalize matrices.
    assert test, "Shift-scale normalization requires a test set."
    assert train0.ncol() and train1.ncol(), "Training samples missing."
    assert test.ncol(), "No test samples found."
    assert normref.ncol(), "No normalization reference samples found."
    assert_rows_aligned(train0, train1, test, normref)

    print("Normalizing with shift-scale, "
          "using a normalization reference with %d samples." % normref.ncol())
    
    # Shift-scale will normalize X to Y.  X is the test data and Y is
    # the training data.
    Y = [None] * train0.nrow()
    for i in range(len(Y)):
        Y[i] = train0._X[i] + train1._X[i]
    DATA_Y = Matrix.InMemoryMatrix(Y)

    # Normalize one sample at a time.
    for j in range(test.ncol()):
        print "Shift-scale normalizing sample %d of %d." % (j+1, test.ncol())
        sys.stdout.flush()
        sample = [x[j] for x in test._X]
        
        X = [None] * normref.nrow()
        for i in range(len(X)):
            X[i] = [sample[i]] + normref._X[i]
        DATA_X = Matrix.InMemoryMatrix(X)
        assert DATA_X.ncol() == normref.ncol() + 1

        DATA_n = shiftscalenorm.normalize(
            DATA_X, DATA_Y, matlab=matlab, binreg_path=binreg_path)
        assert DATA_n.dim() == DATA_X.dim()

        # For debugging:
        #arrayio.tdf.write(DATA_X, open("SAMPLE%02d_nonorm.out"%j, 'w'))
        #arrayio.tdf.write(DATA_n, open("SAMPLE%02d_norm.out"%j, 'w'))

        # Reassign normalized values in place.
        X_n = DATA_n._X
        for i in range(len(X_n)):
            test._X[i][j] = X_n[i][0]

def filter_missing_values(train0, train1, test, normref):
    import math
    
    # Optimization: Assume that the rows are aligned.
    I = []
    for i in range(train0.nrow()):
        x = train0._X[i] + train1._X[i]
        if test:
            x += test._X[i]
        if normref:
            x += normref._X[i]
        x = [x for x in x if math.isnan(x)]
        if not x:
            I.append(i)
        else:
            print "Filtering out nan's found in row %d." % (i+1)
    if len(I) < train0.nrow():
        train0 = train0.matrix(I, None)
        train1 = train1.matrix(I, None)
        if test:
            test = test.matrix(I, None)
        if normref:
            normref = normref.matrix(I, None)
    x = train0, train1, test, normref
    return x
    

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
    
def run_analysis(analysis, options, train0, train1, test, normref):
    import time

    start_time = time.time()

    x = "Analyzing with %d genes and %d metagenes." % (
        analysis.params.genes, analysis.params.metagenes)
    if analysis.params.metagenes == 1:
        x = x.replace("metagenes", "metagene")
    print x

    # Format the parameters and output files for binreg.
    print "Formatting data for binreg analysis."
    file_layout = make_file_layout(options.outpath, analysis.subpath)
    init_paths(file_layout)

    if False and normref:
        # This gives exactly the same results as the standard
        # analysis.  Don't bother running it.  It's slower.
        _run_binreg_normref(
            analysis, options, file_layout, train0, train1, test, normref)
    else:
        _run_binreg_standard(
            analysis, options, file_layout, train0, train1, test)

    # Generate some files for output.
    print "Summarizing results."; sys.stdout.flush()
    summarize_model(file_layout)
    summarize_parameters(analysis.params, file_layout)
    summarize_probabilities(train0, train1, test, file_layout)
    summarize_signature_dataset(file_layout)

    # Make some plots, if desired.
    if not analysis.params.noplots:
        summarize_signature_heatmap(
            options.python, options.arrayplot, options.cluster,
            file_layout, options.libpath)
        summarize_dataset_heatmap(
            options.python, options.arrayplot, options.cluster,
            file_layout, options.libpath)
        summarize_predictions(
            options.povray, analysis.params.label_samples,
            analysis.params.draw_error_bars, file_layout)

    # Make a report for this analysis.
    summarize_report(options.outpath, analysis, start_time)
    sys.stdout.flush()

def _run_binreg_standard(analysis, options, file_layout, train0, train1, test):
    from genomicode import binreg

    print "Running binreg."; sys.stdout.flush()
    write_files_for_binreg(train0, train1, test, file_layout)
    params = analysis.params
    br_params = binreg.BinregParams(
        binreg_version=params.binreg_version,
        cross_validate=int(params.cross_validate),
        make_plots=0,
        num_genes=params.genes, num_metagenes=params.metagenes,
        quantile_normalize=0, shift_scale_normalize=0,
        num_burnin=params.burnin, num_iterations=params.samples,
        num_skips=params.skips,
        credible_interval=params.credible_interval)
    r = binreg.binreg_raw(
        file_layout.BR_EXPRESSION, file_layout.BR_DESCRIPTION,
        1, br_params, matlab=options.matlab,
        binreg_path=options.binreg_path, outpath=file_layout.BINREG)
    print r.read()
    sys.stdout.flush()
    
def _run_binreg_normref(
    analysis, options, file_layout, train0, train1, test, normref):
    import shutil
    from genomicode import Matrix
    from genomicode import binreg
    from genomicode import filelib
    import arrayio

    # test must be specified to shiftscale normalize matrices.
    assert test, "Shift-scale normalization requires a test set."
    assert train0.ncol() and train1.ncol(), "Training samples missing."
    assert test.ncol(), "No test samples found."
    assert normref.ncol(), "No normalization reference samples found."
    assert_rows_aligned(train0, train1, test, normref)

    print("Running binreg "
          "using a normalization reference with %d samples." % normref.ncol())

    for j in range(test.ncol()):
        print "Running binreg on sample %d of %d." % (j+1, test.ncol())
        sys.stdout.flush()

        # Put the test sample in the context of the normalization reference.
        X = [None] * test.nrow()
        for i in range(test.nrow()):
            X[i] = [test._X[i][j]] + normref._X[i]
        x = Matrix.InMemoryMatrix(X, row_names=test._row_names)
        synonyms = { arrayio.ROW_ID : test._synonyms[arrayio.ROW_ID] }
        test_j = Matrix.add_synonyms(x, synonyms)
        
        write_files_for_binreg(train0, train1, test_j, file_layout)
    
        params = analysis.params
        br_params = binreg.BinregParams(
            binreg_version=params.binreg_version,
            cross_validate=int(params.cross_validate),
            make_plots=0,
            num_genes=params.genes, num_metagenes=params.metagenes,
            quantile_normalize=0, shift_scale_normalize=0,
            num_burnin=params.burnin, num_iterations=params.samples,
            num_skips=params.skips,
            credible_interval=params.credible_interval)
        r = binreg.binreg_raw(
            file_layout.BR_EXPRESSION, file_layout.BR_DESCRIPTION,
            1, br_params, matlab=options.matlab,
            binreg_path=options.binreg_path, outpath=file_layout.BINREG)
        br_output = r.read()

        assert os.path.exists(file_layout.BR_VALIDATION), \
               "I could not find the validation file."
        sample_file = "%s.%04d" % (file_layout.BR_VALIDATION, j)
        shutil.move(file_layout.BR_VALIDATION, sample_file)
        
    # Print out last one, so user can see if there are any errors.
    print br_output
    sys.stdout.flush()

    # Combine the validation cases file.
    data = []  # list of (index, type, prob, lower_ci, upper_ci, mgene)
    for j in range(test.ncol()):
        sample_file = "%s.%04d" % (file_layout.BR_VALIDATION, j)
        assert os.path.exists(sample_file)
        iter = filelib.read_cols(sample_file)
        x = iter.next()
        assert len(x) == 6
        data.append(x)
    # Fix the indexes.
    assert data and len(data) == test.ncol()
    first_index = int(data[0][0])
    for j in range(test.ncol()):
        x = data[j]
        x[0] = first_index + j
        data[j] = x
    # Write out the new validationcases file.
    assert not os.path.exists(file_layout.BR_VALIDATION)
    handle = open(file_layout.BR_VALIDATION, 'w')
    for x in data:
        print >>handle, "\t".join(map(str, x))
    handle.close()
        
    

def summarize_model(file_layout):
    from genomicode import filelib
    
    assert filelib.exists_nz(file_layout.BR_COEFFICIENTS), \
           "Cannot find Binreg model.  Binreg failed?"

    handle = open(file_layout.MODEL, 'w')
    print >>handle, "%s\t%s" % ("Name", "Coefficient")
    for line in open(file_layout.BR_COEFFICIENTS):
        x = line.strip().split()
        assert len(x) == 2
        coef, gene = x
        x = gene, coef
        print >>handle, "\t".join(x)
    handle.close()

def summarize_parameters(params, file_layout):
    handle = open(file_layout.PARAMETERS, 'w')
    print >>handle, "NAME\tVALUE"
    print >>handle, "Binreg Version\t%d" % params.binreg_version
    print >>handle, "Genes\t%d" % params.genes
    print >>handle, "Metagenes\t%d" % params.metagenes

    print >>handle, "Train0\t%s" % os.path.split(params.train0_file)[1]
    print >>handle, "Train1\t%s" % os.path.split(params.train1_file)[1]
    x = ""
    if params.test_file:
        x = os.path.split(params.test_file)[1]
    print >>handle, "Test\t%s" % x
    x = ""
    if params.normref_file:
        x = os.path.split(params.normref_file)[1]
    print >>handle, "Normalization Reference File\t%s" % x
    print >>handle, "Log Train0\t%d" % int(params.log_train0)
    print >>handle, "Log Train1\t%d" % int(params.log_train1)
    print >>handle, "Log Test\t%d" % int(params.log_test)
    print >>handle, "Log Normalization Reference\t%d" % int(params.log_normref)
    print >>handle, "Strip AFFX control\t%d" % int(params.strip_affx)
    
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
    print >>handle, "Label Samples\t%d" % int(params.label_samples)
    print >>handle, "Draw Error Bars\t%d" % int(params.draw_error_bars)
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
    assert filelib.exists_nz(file_layout.BR_COEFFICIENTS)
    x = [x[1] for x in filelib.read_cols(file_layout.BR_COEFFICIENTS)]
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
    from genomicode import graphlib

    DATA = arrayio.gct_format.read(file_layout.SIGNATURE_GCT)
    nrow, ncol = DATA.dim()

    # Make the heatmap big enough to read the gene names.  SCALE 1.2
    # is probably the minimum size, but the names are a bit fuzzy.
    SCALE = 1.5
    x = graphlib.find_tall_heatmap_size(
        nrow, ncol, min_box_height=1, min_box_width=1,
        max_box_height=40, max_box_width=40, max_total_height=768*SCALE,
        max_total_width=1024*SCALE)
    xpix, ypix = x
    #print ypix, xpix, nrow, ncol; sys.stdout.flush()

    graphlib.plot_heatmap(
        file_layout.SIGNATURE_GCT, file_layout.SIGNATURE_PNG, xpix, ypix,
        color="bild", show_colorbar=True, show_grid=True,
        gene_center="mean", gene_normalize="var",
        gene_label=True, array_label=True,
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

def summarize_dataset_heatmap(
    python, arrayplot, cluster, file_layout, libpath=[]):
    import arrayio
    from genomicode import graphlib

    DATA = arrayio.gct_format.read(file_layout.DS_SIG)
    nrow, ncol = DATA.dim()

    # 10 won't work for large data sets (e.g. 1000 samples).
    #min_box_width = 10
    min_box_width = 1
    x = graphlib.find_tall_heatmap_size(
        nrow, ncol, min_box_width=min_box_width,
        max_total_height=768*3, max_total_width=1024*3)
    xpix, ypix = x
    #print xpix, ypix, ncol, nrow; sys.stdout.flush()

    graphlib.plot_heatmap(
        file_layout.DS_SIG, file_layout.DS_SIG_PNG, xpix, ypix,
        color="bild", show_colorbar=True, show_grid=True, 
        cluster_genes=True, gene_center="mean", gene_normalize="var",
        array_label=True, 
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

def summarize_predictions(povray, label_samples, draw_error_bars, file_layout):
    # May not plot anything if there are problems with the BinReg
    # output (e.g. all predictions are nan).
    import math
    from genomicode import filelib
    from genomicode import graphlib
    from genomicode import graphconst

    assert filelib.exists_nz(file_layout.PROBABILITIES)

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
        x = float(d.Metagene)
        y = float(d.Probability)
        if math.isnan(x) or math.isnan(y):
            continue

        sample.append(d.Sample)
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
        p = graphconst.SQUARE
        if d.Type in ["train0", "train1"]:
            p = graphconst.CIRCLE

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

    xtick = graphlib.place_ticks(min(X), max(X))
    xtick_label = True
    ytick = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    ytick_label = ["%d%%" % int(x*100) for x in ytick]
    plot_width, plot_height = 1024, 768
    plot_width, plot_height = int(plot_width*1.5), int(plot_height*1.5)
    font = os.path.join(
        os.path.split(graphlib.__file__)[0], "Verdana Bold.ttf")
    onpoint_label = None   # Don't display labels.
    overpoint_label = sample
    if not label_samples:
        overpoint_label = None
    if not draw_error_bars:
        error_bar = None
    points = zip(X, Y)
    graph = graphlib.scatter(
        points, color=color, shape=pch, error_bar=error_bar, point_size=1,
        onpoint_label=onpoint_label, overpoint_label=overpoint_label,
        ylim=(0, 1.02), 
        xtick=xtick, xtick_label=xtick_label,
        ytick=ytick, ytick_label=ytick_label, tick_size=1, 
        xlabel="Metagene Score", ylabel="Probability",
        label_size=1, width=plot_width, height=plot_height,
        font=font)
    x = graph.write(file_layout.PREDICTIONS_PNG, povray_bin=povray)
    print x; sys.stdout.flush()
    assert filelib.exists_nz(file_layout.PREDICTIONS_PNG), \
           "Failed to plot predictions."
    # povray -D -J +Opredictions.png -H768 -W1024 +A0.5 predictions.pov

def summarize_report(outpath, analysis, start_time):
    import time
    from genomicode import parselib
    from genomicode import htmllib

    params = analysis.params
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")

    lines = []
    w = lines.append
    w("<HTML>")
    w(htmllib.HEAD(htmllib.TITLE("CreateSignatures Report")))
    w("<BODY>")
    w(htmllib.CENTER(htmllib.H1(htmllib.EM("CreateSignatures") + " Report")))

    # Describe the analysis in human language.
    w(htmllib.H3("I.  Analysis"))

    items = []

    # Files.
    t0 = os.path.split(params.train0_file)[1]
    t1 = os.path.split(params.train1_file)[1]
    te = None
    nr = None
    if params.test_file:
        te = os.path.split(params.test_file)[1]
    if params.normref_file:
        nr = os.path.split(params.normref_file)[1]
    x = "I ran a signature analysis using a training set of %s and %s." % (
        highlight(t0), highlight(t1))
    if te:
        x += "  I generated predictions on %s." % highlight(te)
    items.append(x)

    # Algorithm.
    genes = highlight("%d genes" % params.genes)
    metagenes = highlight("%d metagenes" % params.metagenes)
    if params.metagenes == 1:
        metagenes = metagenes.replace("genes", "gene")
    x = "I used the %s algorithm with %s and %s." % (
        highlight("BinReg %d" % params.binreg_version), genes, metagenes)
    items.append(x)
        
    # Logging the data.
    logged = []
    if params.log_train0:
        logged.append(t0)
    if params.log_train1:
        logged.append(t1)
    if params.log_test:
        assert te
        logged.append(te)
    if params.log_normref:
        assert nr
        logged.append(nr)
    # Make sure there are no duplicates in logged.
    logged = [logged[i] for i in range(len(logged))
              if logged[i] not in logged[:i]]
    logged_str = None
    if len(logged) == 1:
        logged_str = logged[0]
    elif len(logged) == 2:
        logged_str = "%s and %s" % tuple(logged)
    elif len(logged) == 3:
        logged_str = "%s, %s, and %s" % tuple(logged)
    elif len(logged) == 4:
        logged_str = "%s, %s, %s, and %s" % tuple(logged)
    elif len(logged):
        raise AssertionError, "Too many logged files."
    if logged_str:
        x = "I logged the expression values in %s." % logged_str
        items.append(x)

    # Normalization.
    norm = []
    if params.quantile:
        norm.append("quantile")
    if params.shiftscale:
        norm.append("Shift-Scale")
    if params.dwd:
        norm.append("DWD")
    if params.dwd_bild:
        norm.append("DWD (Bild method)")
    norm_str = None
    if len(norm) == 1:
        norm_str = highlight(norm[0])
    elif len(norm) == 2:
        norm_str = "%s and %s" % tuple(map(highlight, norm))
    elif len(norm) == 3:
        norm_str = "%s, %s, and %s" % tuple(map(highlight, norm))
    elif len(norm) == 4:
        norm_str = "%s, %s, %s, and %s" % tuple(map(highlight, norm))
    if norm_str:
        x = "I applied %s normalization." % norm_str
    else:
        x = "I did not normalize the data."
    if params.normref_file:
        if params.shiftscale:
            nrf = os.path.split(params.normref_file)[1]
            x += "  For shift-scale, I used %s as a normalization reference."\
                 % nrf
    items.append(x)
        
    # Strip Affymetrix control probes.
    if params.strip_affx:
        x = "As requested, I stripped the Affymetrix control probes."
        items.append(x)

    # MCMC
    x = ("For the statistical (Markov chain Monte Carlo) simulation, I "
         "discarded %s samples for the burn-in and then collected %s "
         "samples for the model." % (
             parselib.pretty_int(params.burnin),
             parselib.pretty_int(params.samples)))
    if params.skips > 1:
        x += "  During the sampling, I recorded every %s sample." % (
            parselib.pretty_ordinal(params.skips))
    items.append(x)

    items = [htmllib.LI() + x for x in items]
    items_str = "\n".join(items)
    w(htmllib.UL(items_str))

    # Put together a table with the results.
    w(htmllib.P())
    w(htmllib.H3("II.  Results"))

    rows = []

    file_layout = make_file_layout(outpath, analysis.subpath)
    x = read_matrices(
        params.train0_file, params.train1_file, params.test_file,
        params.normref_file)
    train0, train1, test, normref = x
    x = _make_signature_figure(file_layout.SIGNATURE_PNG, train0, train1)
    figure, legend = x
    figure_col1 = htmllib.TD(figure)
    legend_col1 = htmllib.TD(legend, valign="TOP")

    x = _make_predictions_figure(
        file_layout.PREDICTIONS_PNG, file_layout.PROBABILITIES, params)
    figure, legend = x
    figure_col2 = htmllib.TD(figure)
    legend_col2 = htmllib.TD(legend, valign="TOP")

    x = htmllib.TR(figure_col1+"\n"+figure_col2)
    rows.append(x)
    x = htmllib.TR(legend_col1+"\n"+legend_col2)
    rows.append(x)
    rows_str = "\n".join(rows)
    w(htmllib.TABLE(rows_str, border=0, cellspacing=10))

    # Write the footer.
    w(htmllib.P())
    w(htmllib.HR())
    end_time = time.time()
    w(_make_footer(start_time, end_time))
    w("</BODY>")
    w("</HTML>")

    x = "\n".join(lines) + "\n"
    outfile = file_layout.REPORT
    open(outfile, 'w').write(x)

def summarize_all_report(outpath, analyses, start_time):
    import time
    from genomicode import htmllib

    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def yesno(x):
        if x:
            return "Yes"
        return "No"

    lines = []
    w = lines.append
    w("<HTML>")
    w(htmllib.HEAD(htmllib.TITLE("CreateSignatures Report")))
    w("<BODY>")
    w(htmllib.CENTER(htmllib.H1(htmllib.EM("CreateSignatures") + " Report")))

    # Put together a table with the summary.
    w(htmllib.P())
    w(htmllib.H3("I ran %d analyses:" % len(analyses)))
    rows = []

    # Make a header.
    cols = [
        htmllib.TD("#"),
        htmllib.TH("Version"),
        htmllib.TH("Genes"),
        htmllib.TH("Metagenes"),
        htmllib.TH("Quantile"),
        htmllib.TH("Shift-Scale"),
        htmllib.TH("DWD"),
        htmllib.TH("DWD (Bild)"),
        htmllib.TH("Train 0", align="LEFT"),
        htmllib.TH("Train 1", align="LEFT"),
        htmllib.TH("Test", align="LEFT"),
        ]
    x = htmllib.TR("\n".join(cols))
    rows.append(x)

    for i, analysis in enumerate(analyses):
        params = analysis.params
        cols = []
        x = htmllib.A("%d." % (i+1), "#analysis%02d" % (i+1))
        cols.append(x)
        cols.append(params.binreg_version)
        cols.append(params.genes)
        cols.append(params.metagenes)
        cols.append(yesno(params.quantile))
        cols.append(yesno(params.shiftscale))
        cols.append(yesno(params.dwd))
        cols.append(yesno(params.dwd_bild))
        cols.append(os.path.split(params.train0_file)[1])
        cols.append(os.path.split(params.train1_file)[1])
        x = "None"
        if params.test_file:
            x = os.path.split(params.test_file)[1]
        cols.append(x)

        for i in range(len(cols)):
            if i >= 0  and i < 8:
                cols[i] = htmllib.TD(cols[i], align="MIDDLE")
            else:
                cols[i] = htmllib.TD(cols[i])
        x = htmllib.TR("\n".join(cols))
        rows.append(x)
    rows_str = "\n".join(rows)
    w(htmllib.TABLE(rows_str, border=1, cellpadding=3, cellspacing=0))


    # Put together a table with each of the results.
    w("<P>\n")
    rows = []

    for i, analysis in enumerate(analyses):
        params = analysis.params
        file_layout = make_file_layout(outpath, analysis.subpath)
        
        if i > 0:
            # Make a divider between the rows.
            x = htmllib.TR(htmllib.TD(htmllib.HR(), colspan=2))
            rows.append(x)

        # Make the heading for this set of parameters.
        x = htmllib.TR(
            htmllib.TD(
            htmllib.A(
                htmllib.H3("%d Genes, %d Metagenes" % (
                        params.genes, params.metagenes))
                    , name="analysis%02d" % (i+1)),
                colspan=2, align="CENTER"))
        rows.append(x)

        x = read_matrices(
            params.train0_file, params.train1_file, params.test_file,
            params.normref_file)
        train0, train1, test, normref = x

        x = _make_signature_figure(
            file_layout.GLOBAL_SIGNATURE_PNG, train0, train1)
        figure, legend = x
        figure_col1 = htmllib.TD(figure)
        legend_col1 = htmllib.TD(legend, valign="TOP")

        x = _make_predictions_figure(
            file_layout.GLOBAL_PREDICTIONS_PNG,
            file_layout.GLOBAL_PROBABILITIES, params)
        figure, legend = x
        figure_col2 = htmllib.TD(figure)
        legend_col2 = htmllib.TD(legend, valign="TOP")

        x = htmllib.TR(figure_col1+figure_col2)
        rows.append(x)
        x = htmllib.TR(legend_col1+legend_col2)
        rows.append(x)
        
    rows_str = "\n".join(rows)
    w(htmllib.TABLE(rows_str, border=0, cellspacing=10))

    # Write the footer.
    w(htmllib.P())
    w(htmllib.HR())
    end_time = time.time()
    w(_make_footer(start_time, end_time))
    w("</BODY>")
    w("</HTML>")

    x = "\n".join(lines) + "\n"
    file_layout = make_file_layout(outpath, analyses[0].subpath)
    outfile = file_layout.GLOBAL_REPORT
    open(outfile, 'w').write(x)

def _make_signature_figure(png_file, train0, train1):
    from genomicode import filelib
    from genomicode import htmllib
    
    x = "No plot generated."
    if filelib.exists_nz(png_file):
        sig_file = os.path.split(png_file)[1]
        x = htmllib.A(htmllib.IMG(height=480, src=sig_file), href=sig_file)
    figure = x

    legend = (
        htmllib.B("Figure 1: Signature Heatmap. ") +
        "In this heatmap, each row represents a gene in the signature. "
        "The first %d columns are the samples from the %s "
        "data set, and the remaining %d columns are the samples from the "
        "%s data set. "
        "Warm colors indicate high expression of the gene, and cool "
        "colors indicate low expression." % (
            train0.ncol(), htmllib.EM("train0"),
            train1.ncol(), htmllib.EM("train1"),
            ))
    return figure, legend
    
def _make_predictions_figure(png_file, probabilities_file, params):
    from genomicode import filelib
    from genomicode import htmllib
    
    x = "No plot generated."
    if filelib.exists_nz(png_file):
        pred_file = os.path.split(png_file)[1]
        x = htmllib.A(
            htmllib.IMG(height=480, src=pred_file), href=pred_file)
    figure = x
    
    x = (
        htmllib.B("Figure 2: Predictions. ") +
        "This scatter plot shows the predictions from the signature "
        "for each sample. "
        "On the Y-axis, high probabilities indicate that the gene "
        "expression profile of the sample better resembles the train 1 "
        "class, while low probabilities indicate a closer resemblance "
        "to train 0.  "
        "The blue and red circles are the predictions (from "
        "leave-one-out cross-validation) on the train 0 and train 1 "
        "samples, respectively. ")
    if params.test_file:
        x += "The black squares are the predictions on the test samples. "
    if params.draw_error_bars:
        x += "The error bars show the %d%% credible interval. " % \
             params.credible_interval
    x += (
        "The X-axis, the Metagene Score, is the magnitude of the "
        "sample on the first principal component. "
        "This is used to help separate the samples across the plot.")
    prob_file = os.path.split(probabilities_file)[1]
    x += htmllib.P() + (
        "The raw values from this plot are available as a "
        'tab-delimited text file: %s.' %
        htmllib.A(prob_file, href=prob_file))
    legend = x

    return figure, legend

def _make_footer(start_time, end_time):
    from genomicode import parselib
    from genomicode import htmllib
    
    time_str = parselib.pretty_date(start_time)
    run_time = pretty_runtime(start_time, end_time)
    hostname = pretty_hostname()
    footer = htmllib.EM(
        "This analysis was run on %s on %s.  It took %s to complete.\n" %
        (time_str, hostname, run_time))
    return footer

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

def pretty_runtime(start_time, end_time):
    from genomicode import parselib
    
    x = end_time-start_time
    fracs = x - int(x)
    fracs = int(fracs * 1000)
    x = int(x)
    num_hours = x / 3600
    x = x % 3600
    num_secs = x % 60
    num_mins = x / 60
    if num_hours == 0 and num_mins == 0:
        run_time = "%ss" % parselib.pretty_int(num_secs)
    elif num_hours == 0:
        run_time = "%02d:%02d.%03d" % (num_mins, num_secs, fracs)
    else:
        run_time = "%02d:%02d:%02d.%03d" % (
            num_hours, num_mins, num_secs, fracs)
    return run_time

def pretty_hostname():
    import subprocess
    
    cmd = "hostname"
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    wh, r = p.stdin, p.stdout
    wh.close()
    hostname = r.read().strip()
    assert hostname, "I could not get the hostname."
    return hostname
    
def main():
    import time

    # Print the starting message.
    start_time = time.time()
    x1 = pretty_hostname()
    x2 = time.strftime("%A, %d %B %Y, %I:%M %p", time.localtime(start_time))
    print "Running pybinreg.py analysis on %s on %s." % (x1, x2)

    print "Reading input parameters."; sys.stdout.flush()
    options, param_list = init_params()
    analyses = make_analyses(param_list)
    assert analyses

    # To DEBUG report.
    #if len(analyses) > 1:
    #    summarize_all_report(options.outpath, analyses, start_time)
    #sys.exit(0)

    # import modules after init_params, if possible.  GenePattern
    # messes up the PATH variable, so many imports won't work until
    # init_params fixes it.
    from genomicode import archive
    from genomicode import parselib
    import arrayio

    print "Setting up a folder for the output files."; sys.stdout.flush()
    # Should set up a file layout for one of the analyses subpaths.
    # If I use None, then will create extra directories
    # (e.g. "binreg") that may not be necessary for a global analysis.
    file_layout = make_file_layout(options.outpath, analyses[0].subpath)
    init_paths(file_layout)

    # Load the files and standardize them to GCT format.
    print "Reading gene expression matrices."; sys.stdout.flush()
    # Assume all analyses have the same train0, train1, test, and
    # normref files.
    params = analyses[0].params
    x = read_matrices_all(
        params.train0_file, params.train1_file, params.test_file,
        params.normref_file)
    t0, t1, te, nr, train0, train1, test, normref = x
    print "train0 file has %s genes and %s samples." % (
        tuple(map(parselib.pretty_int, t0.dim())))
    print "train1 file has %s genes and %s samples." % (
        tuple(map(parselib.pretty_int, t1.dim())))
    if te:
        print "test file has %s genes and %s samples." % (
            tuple(map(parselib.pretty_int, te.dim())))
    if nr:
        print "normalization reference file has %s genes and %s samples." % (
            tuple(map(parselib.pretty_int, nr.dim())))
    print "Merged file has %s genes." % (parselib.pretty_int(train0.nrow()))
    train0 = arrayio.convert(train0, to_format=arrayio.gct_format)
    train1 = arrayio.convert(train1, to_format=arrayio.gct_format)
    if test:
        test = arrayio.convert(test, to_format=arrayio.gct_format)
    if normref:
        normref = arrayio.convert(normref, to_format=arrayio.gct_format)
    write_dataset(file_layout.DS_ORIG, train0, train1, test)
    if normref:
        arrayio.gct_format.write(normref, open(file_layout.NR_ORIG, 'w'))
    sys.stdout.flush()

    # Remove the Affymetrix control probes.  Do this before any
    # normalization.
    if params.strip_affx:
        x = strip_affx_control_probes(train0, train1, test, normref)
        train0, train1, test, normref = x

    # Log normalize each of the files, if necessary.  Will log in place.
    logged = log_matrices(
        train0, train1, test, normref, params.log_train0, params.log_train1,
        params.log_test, params.log_normref)
    if logged:
        write_dataset(file_layout.DS_LOG, train0, train1, test)

    # Apply the appropriate normalizations.
    if params.quantile:
        quantile_normalize_matrices(
            train0, train1, test, normref, matlab=options.matlab,
            binreg_path=options.binreg_path)
        write_dataset(file_layout.DS_QNORM, train0, train1, test)
        if normref:
            arrayio.gct_format.write(normref, open(file_layout.NR_QNORM, 'w'))

    if params.dwd:
        dwd_normalize_matrices(
            train0, train1, test, matlab=options.matlab,
            dwd_path=options.dwd_path)
        write_dataset(file_layout.DS_DWD, train0, train1, test)

    if params.dwd_bild:
        dwd_normalize_matrices(
            train0, train1, test, version="bild", matlab=options.matlab,
            dwd_path=options.dwd_path)
        write_dataset(file_layout.DS_DWD_BILD, train0, train1, test)

    if params.shiftscale:
        shiftscale_normalize_matrices(
            train0, train1, test, normref, matlab=options.matlab,
            binreg_path=options.binreg_path)
        write_dataset(file_layout.DS_SS, train0, train1, test)

    # Make sure there are no NA's, or BinReg will not do the
    # predictions correctly.
    x = filter_missing_values(train0, train1, test, normref)
    train0, train1, test, normref = x

    print "Writing processed dataset."; sys.stdout.flush()
    write_dataset(file_layout.DS_FINAL, train0, train1, test)
    if normref:
        arrayio.gct_format.write(normref, open(file_layout.NR_FINAL, 'w'))

    # Run each of the binreg analyses.
    for i, analysis in enumerate(analyses):
        if len(analyses) > 1:
            print "Running analysis %d of %d." % (i+1, len(analyses))
            sys.stdout.flush()
        run_analysis(analysis, options, train0, train1, test, normref)

    # If there were multiple analyses specified, make a copy of
    # some files for convenience.  Need to copy it rather than
    # symlink, in case we archive it later.
    if len(analyses) > 1:
        for analysis in analyses:
            file_layout = make_file_layout(options.outpath, analysis.subpath)

            files_to_copy = [
                (file_layout.PROBABILITIES, file_layout.GLOBAL_PROBABILITIES),
                (file_layout.PREDICTIONS_PNG,
                 file_layout.GLOBAL_PREDICTIONS_PNG),
                (file_layout.SIGNATURE_PNG, file_layout.GLOBAL_SIGNATURE_PNG),
                ]
            for src, dst in files_to_copy:
                assert os.path.exists(src)
                os.system("cp -p '%s' '%s'" % (src, dst))

            if analysis.params.archive:
                archive.zip_path(file_layout.ANALYSIS)
        summarize_all_report(options.outpath, analyses, start_time)

    # Assume all archive options are the same.
    params = analyses[0].params
    if params.archive:
        print "Compressing results."; sys.stdout.flush()
        file_layout = make_file_layout(options.outpath, analyses[0].subpath)
        if os.path.exists(file_layout.BINREG):
            archive.zip_path(file_layout.BINREG)
        if os.path.exists(file_layout.ATTIC):
            archive.zip_path(file_layout.ATTIC)
        if os.path.exists(file_layout.GLOBAL_ATTIC):
            archive.zip_path(file_layout.GLOBAL_ATTIC)

    #make_FILES(files)
    
    x = time.strftime("%A, %d %B %Y, %I:%M %p", time.localtime())
    print "Finished pybinreg.py analysis on %s." % x

if __name__ == '__main__':
    main()
    #import profile; profile.run("main()")
