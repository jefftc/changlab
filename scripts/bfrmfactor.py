#!/usr/bin/env python

import os, sys

# make_file_layout
# init_paths
# 
# write_dataset
# write_bfrm_dataset
# write_bfrm_bin
# write_ids
# write_factor_ids
# write_h
# write_evol
#
# run_bfrm
# log_matrix
# filter_dataset
# create_control_vars
#
# summarize_factor_scores
# summarize_gene_factor_probs
# summarize_factor_geneset
#
# _read_model

def make_file_layout(outpath):
    from genomicode import filelayout

    outpath = outpath or "."
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath = os.path.realpath(outpath)

    Path, File = filelayout.Path, filelayout.File
    
    # Have only set set of these files for the whole analysis.
    FILES = [
        Path.ATTIC("attic",
            File.DATASET_ORIG("dataset_orig.gct"),
            ),
        File.DATASET("dataset.gct"),
        File.FACTOR_SCORES("factors.pcl"),
        File.FACTOR_CDT("factors.cdt"),   # created by cluster
        File.FACTOR_ATR("factors.atr"),   # created by cluster
        File.FACTOR_GTR("factors.gtr"),   # created by cluster
        File.FACTOR_SCORES_PNG("factors.png"),
        File.FACTOR_PROBS("p_gene_factor.txt"),
        File.FACTOR_PROBS_PNG("p_gene_factor.png"),
        File.FACTOR_PROBS_ALL("p_geneall_factor.txt"),
        File.FACTOR_GENESET("genesets.gmt"),

        Path.BFRM("bfrm_model",
            File.BFRM_BIN("bfrm64"),
            File.BFRM_PARAMETERS("parameters.txt"),
            File.BFRM_DATASET("dataset.txt"),
            File.BFRM_GENE_IDS("geneids.txt"),
            File.BFRM_PROBE_IDS("probeids.txt"),
            File.BFRM_SAMPLE_IDS("sids.txt"),
            File.BFRM_FACTOR_IDS("factorids.txt"),
            File.BFRM_H("H.txt"),
            File.BFRM_EVOL("evol_i.txt"),
                  
            File.BFRM_MA("mA.txt"),
            ),
        ]
    
    file_layout = Path.OUTPATH(outpath, *FILES)
    return file_layout

def init_paths(file_layout):
    from genomicode import filelayout

    for x in filelayout.walk(file_layout):
        dirpath, dirnames, filenames = x
        if os.path.exists(dirpath):
            continue
        os.mkdir(dirpath)

def write_dataset(filename, MATRIX):
    import arrayio
    arrayio.gct_format.write(MATRIX, open(filename, 'w'))

def write_bfrm_dataset(filename, DATA):
    data = DATA.value()
    handle = open(filename, 'w')
    for x in data:
        print >>handle, "\t".join(map(str, x))
    handle.close()

def write_bfrm_bin(filename, bfrm_bin):
    import stat
    from genomicode import bfrm
    
    bfrm_bin = bfrm.find_bfrm_bin(bfrm_bin)
    assert bfrm_bin and os.path.exists(bfrm_bin)
    x = open(bfrm_bin, 'rb').read()
    open(filename, 'wb').write(x)
    os.chmod(filename, stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)

def write_ids(geneid_file, probeid_file, sampleid_file, DATA):
    import arrayio
    from genomicode import bfrm
    
    # Gene ID file contains the row IDS from the matrix, with a header.
    row_names = DATA.row_names()
    handle = open(geneid_file, 'w')
    print >>handle, "\t".join(row_names)
    for i in range(DATA.nrow()):
        x = [DATA.row_names(x)[i] for x in row_names]
        print >>handle, "\t".join(x)
    handle.close()

    # Probe ID file contains the Affymetrix probeset IDs, no header.
    name = bfrm.get_affy_row_name(DATA)
    probeset_ids = DATA.row_names(name)
    x = ["%s\n" % x for x in probeset_ids]
    open(probeid_file, 'w').writelines(x)

    # Sample ID file contains the samples, with no header.
    sample_ids = DATA.col_names(arrayio.COL_ID)
    x = ["%s\n" % x for x in sample_ids]
    open(sampleid_file, 'w').writelines(x)

def write_factor_ids(filename, file_layout):
    from genomicode import parselib
    
    model = _read_model(file_layout)
    num_factors = len(model["FACTOR_O"])

    if not num_factors:
        # If no factors, then just make an empty file.
        open(filename, 'w')
        return

    # Create some names for the factors.
    x = parselib.pretty_range(1, num_factors+1)
    factor_ids = ["FACTOR%s" % x for x in x]
    
    x = ["%s\n" % x for x in factor_ids]
    open(filename, 'w').writelines(x)

def write_h(filename, CONTROL):
    # sample x design factor matrix
    # 1st column = 1, intercept
    # 2-n column = control vectors
    assert CONTROL
    
    handle = open(filename, 'w')
    for i in range(len(CONTROL[0])):
        x = [1] + [x[i] for x in CONTROL]
        print >>handle, "\t".join(map(str, x))
    handle.close()

def write_evol(filename, DATA, nucleus_file):
    # EVOL contains the 1-based indexes of the genes.
    assert os.path.exists(nucleus_file)
    x = open(nucleus_file).read()
    names = x.strip().split()
    
    x = DATA._index(row=names)
    I_row, I_col = x
    assert I_row is not None
    assert I_row, "I could not find any of the nucleus genes."

    print "I matched %d genes for the evolutionary search." % len(I_row)

    # BFRM needs 1-based indexes.
    indexes = [i+1 for i in I_row]
    handle = open(filename, 'w')
    for i in indexes:
        print >>handle, i
    handle.close()

def run_bfrm(
    file_layout, bfrm_bin, num_control_vars, start_factors, nucleus_file,
    max_factors, max_genes):
    import random
    import time
    import subprocess
    import arrayio
    from genomicode import bfrm

    DATA = arrayio.read(file_layout.DATASET)

    # Configure the parameters file.
    params = bfrm.get_default_parameters(bfrm_bin=bfrm_bin)

    # Copy the BFRM binary.
    write_bfrm_bin(file_layout.BFRM_BIN, bfrm_bin)

    # Write the geneids.txt and sampleids.txt files.
    write_ids(
        file_layout.BFRM_GENE_IDS, file_layout.BFRM_PROBE_IDS,
        file_layout.BFRM_SAMPLE_IDS, DATA)
    
    # Write the data to a file.
    write_bfrm_dataset(file_layout.BFRM_DATASET, DATA)
    dataset = os.path.split(file_layout.BFRM_DATASET)[1]
    params = bfrm.set_params_dataset(
        params, DATA.ncol(), DATA.nrow(), dataset)

    # Set the starting number of latent factors.
    if start_factors:
        x = "Initiating analysis with %d factors." % start_factors
        if start_factors == 1:
            x = x.replace("factors", "factor")
        print x
        assert start_factors < DATA.nrow()
        params = bfrm.set_params_latent_factors(params, start_factors)

    # If control variables are requested, then make an H file.
    if num_control_vars:
        print "Using %d control variables." % num_control_vars
        # Should create control variables based on original data set.
        DATA_orig = arrayio.read(file_layout.DATASET_ORIG)
        CONTROL = create_control_vars(DATA_orig, num_control_vars)
        assert len(CONTROL) == num_control_vars
        assert len(CONTROL[0]) == DATA_orig.ncol()
        write_h(file_layout.BFRM_H, CONTROL)
        h = os.path.split(file_layout.BFRM_H)[1]
        params = bfrm.set_params_control(params, len(CONTROL), h)

    # If the nucleus was provided, then set up an evolutionary search.
    if nucleus_file:
        print "Enabling evolutionary search."
        write_evol(file_layout.BFRM_EVOL, DATA, nucleus_file)
        params = bfrm.set_params_evol(
            params, file_layout.BFRM_EVOL, max_factors, max_genes)

    # Write out the parameters file.
    open(file_layout.BFRM_PARAMETERS, 'w').write(params)

    # Execute BFRM.
    print "Running BFRM."
    sys.stdout.flush()

    # If I run many factor analyses at the same time, they'll generate
    # the same results.  Thus, delay a little bit so as to shuffle
    # around the random number generator.
    NUM_SECONDS = 10   # maximum time to delay
    random.seed(os.getpid()*time.time())
    time.sleep(random.random()*NUM_SECONDS)

    cwd = os.getcwd()
    try:
        os.chdir(file_layout.BFRM)
        cmd = "%s %s" % (file_layout.BFRM_BIN, file_layout.BFRM_PARAMETERS)
        #print cmd
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        for line in r:
            sys.stdout.write(line)
            sys.stdout.flush()
    finally:
        os.chdir(cwd)
    assert os.path.exists(file_layout.BFRM_MA)

    write_factor_ids(file_layout.BFRM_FACTOR_IDS, file_layout)

    # Print out the number of factors found.
    model = _read_model(file_layout)
    num_genes = len(model["GENE_O"])
    num_factors = len(model["FACTOR_O"])
    if not num_factors:
        print "Found %d factors." % (num_factors)
    else:
        x = "Found %d genes and %d factors." % (num_genes, num_factors)
        if num_factors == 1:
            x = x.replace("factors", "factor")
        print x

def log_matrix(MATRIX):
    # Log the matrix if necessary.  Will log in place.  Return a
    # boolean indicating whether anything was logged.
    from genomicode import jmath
    from genomicode import binreg

    if binreg.is_logged_array_data(MATRIX):
        return False
    print "I will log the matrix."
    MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)
    return True

def filter_dataset(MATRIX, filter_mean, filter_var):
    from genomicode import pcalib

    if not filter_mean and not filter_var:
        return MATRIX

    # Calculate the number of genes to keep.
    num_genes = MATRIX.nrow()
    num_genes_mean = num_genes_var = None
    if filter_mean:
        num_genes_mean = int(round(num_genes * (1.0-filter_mean)))
    if filter_var:
        num_genes_var = int(round(num_genes * (1.0-filter_var)))
    I = pcalib.select_genes_mv(MATRIX._X, num_genes_mean, num_genes_var)
    MATRIX_f = MATRIX.matrix(I, None)
    return MATRIX_f

def create_control_vars(MATRIX, num_control_vars):
    import numpy
    from genomicode import jmath
    from genomicode import bfrm
    
    # Look for the annotations that resemble affymetrix probe set IDs.
    affx_name = bfrm.get_affy_row_name(MATRIX)

    # Select the affymetrix control variables.
    ids = MATRIX.row_names(affx_name)
    I = [i for (i, id) in enumerate(ids) if bfrm.is_affx(id)]
    assert I
    AFFX = MATRIX.matrix(I, None)
    max_control_vars = min(AFFX.nrow(), AFFX.ncol())
    assert num_control_vars<= AFFX.nrow() and num_control_vars <= AFFX.ncol(),\
           "Too many control variables.  Maximum is %d." % (
        num_control_vars, max_control_vars)

    # Calculate the SVD of the control probes.
    X = AFFX._X
    # Subtract the means from each gene.
    for i in range(len(X)):
        x = X[i]
        m = jmath.mean(x)
        x = [x-m for x in x]
        X[i] = x
    # Calculate the SVD.
    U, s, V = numpy.linalg.svd(X, full_matrices=False)
    # Each row of V is a control factor.
    CONTROL = V.tolist()[:num_control_vars]
    assert len(CONTROL) == num_control_vars
    assert len(CONTROL[0]) == AFFX.ncol()
    
    return CONTROL

def summarize_factor_scores(
    file_layout, factor_cutoff, python, arrayplot, cluster, libpath):
    import arrayio
    from genomicode import Matrix
    from genomicode import plotlib

    DATA = arrayio.read(file_layout.DATASET)
    model = _read_model(file_layout, factor_cutoff)

    F = model["F"]
    # If there were no factors, then don't generate any files.
    if not F.nrow():
        print "Not generating factor scores file.  No factors detected."
        return
    assert F.ncol() == DATA.ncol()

    # Read the factor names.
    x = [x.strip() for x in open(file_layout.BFRM_FACTOR_IDS)]
    factor_names = x
    assert len(factor_names) == F.nrow()
    # The factor names are in the same order as the data files.  Sort
    # them so they'll be in the same order as the clean model.
    factor_names = [factor_names[i] for i in model["FACTOR_O"]]

    SAMPLE_NAME = arrayio.tdf.SAMPLE_NAME
    row_names = {}
    col_names = {}
    row_names["xID"] = factor_names
    col_names[SAMPLE_NAME] = DATA.col_names(SAMPLE_NAME)
    M = Matrix.InMemoryMatrix(F._X, row_names, col_names)
    arrayio.pcl_format.write(M, file_layout.FACTOR_SCORES)

    # Make the heatmap.
    x = plotlib.find_wide_heatmap_size(
        M.nrow(), M.ncol(), min_box_height=10, min_box_width=10,
        max_total_height=768, max_total_width=1024)
    xpix, ypix = x
    ypix = min(ypix, xpix*4)
    x = plotlib.plot_heatmap(
        file_layout.FACTOR_SCORES, file_layout.FACTOR_SCORES_PNG,
        xpix, ypix, gene_label=True, cluster_genes=True,
        gene_center="mean", gene_normalize="var", array_label=True,
        cluster_arrays=True,
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

    # Clean up some of the cluster files.
    files = [
        file_layout.FACTOR_CDT, file_layout.FACTOR_ATR, file_layout.FACTOR_GTR
        ]
    for filename in files:
        if not os.path.exists(filename):
            continue
        src = filename
        x = os.path.split(filename)[1]
        dst = os.path.join(file_layout.ATTIC, x)
        os.rename(src, dst)

def summarize_gene_factor_probs(
    file_layout, factor_cutoff, python, arrayplot, cluster, libpath):
    import arrayio
    from genomicode import Matrix
    from genomicode import plotlib

    model = _read_model(file_layout, factor_cutoff)
    PostPib = model["PostPib"]
    ExternalProb = model.get("ExternalProb")

    # If there were no factors, then don't generate any files.
    if not PostPib.ncol():
        print "Not generating factor probabilities file.  No factors detected."
        return

    # Pull out the gene names.
    DATA = arrayio.read(file_layout.DATASET)
    DATA_m = DATA.matrix(model["VariablesIn"], None)

    # Pull out the factor names.
    assert os.path.exists(file_layout.FACTOR_SCORES)
    D_scores = arrayio.read(file_layout.FACTOR_SCORES)
    factor_names = D_scores.row_names(arrayio.ROW_ID)
    assert len(factor_names) == PostPib.ncol()

    # Write the probabilities for the genes in the model.
    SAMPLE_NAME = arrayio.tdf.SAMPLE_NAME
    row_names = {}
    col_names = {}
    row_order = DATA_m.row_names()
    for x in row_order:
        row_names[x] = DATA_m.row_names(x)
    col_names[SAMPLE_NAME] = factor_names
    M = Matrix.InMemoryMatrix(PostPib._X, row_names, col_names, row_order)
    arrayio.tab_delimited_format.write(M, file_layout.FACTOR_PROBS)

    # Make heatmap of the factor probs.
    #x = plotlib.find_tall_heatmap_size(
    #    M.nrow(), M.ncol(), min_box_width=10, max_total_height=1000,
    #    max_total_width=1000)
    xpix, ypix = 20, 20
    x = plotlib.plot_heatmap(
        file_layout.FACTOR_PROBS, file_layout.FACTOR_PROBS_PNG, xpix, ypix,
        color="red", array_label=True, scale=-0.5, gain=2.0,
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)
    
    # If exists, write the probabilities for all genes in the data set.
    if not ExternalProb:
        return
    row_names = {}
    col_names = {}
    row_order = DATA.row_names()
    for x in row_order:
        row_names[x] = DATA.row_names(x)
    col_names[SAMPLE_NAME] = factor_names
    M = Matrix.InMemoryMatrix(ExternalProb._X, row_names, col_names, row_order)
    arrayio.tab_delimited_format.write(M, file_layout.FACTOR_PROBS_ALL)

def summarize_factor_geneset(file_layout, factor_cutoff):
    import arrayio

    # FACTOR_PROBS should be made first.  It may not exist if there
    # are no factors.
    if not os.path.exists(file_layout.FACTOR_PROBS):
        return
    M = arrayio.read(file_layout.FACTOR_PROBS)

    # Choose the right IDs for the gene.  Should be in GCT format.
    name = "DESCRIPTION"
    if name not in M.row_names():
        name = arrayio.ROW_ID
    gene_names = M.row_names(name)

    factor2genes = {}
    factor_names = M.col_names(arrayio.tdf.SAMPLE_NAME)
    assert len(factor_names) == M.ncol()
    for i in range(M.ncol()):
        probs = M.value(None, i)
        assert len(probs) == len(gene_names)
        genes = []
        for p, n in zip(probs, gene_names):
            if p >= factor_cutoff:
                genes.append(n)
        genes = [x for x in genes if x]
        genes = [x for x in genes if x != "---"]
        genes.sort()
        factor2genes[factor_names[i]] = genes

    # Write out the factors in GMT format.
    handle = open(file_layout.FACTOR_GENESET, 'w')
    for fn in factor_names:
        genes = factor2genes[fn]
        x = [fn, "na"] + genes
        print >>handle, "\t".join(x)
    handle.close()

def _read_model(file_layout, factor_cutoff=None):
    from genomicode import bfrm
    
    param_file = os.path.split(file_layout.BFRM_PARAMETERS)[1]
    model = bfrm.read_clean_model(
        file_layout.BFRM, param_file=param_file, factor_cutoff=factor_cutoff)
    return model

def main():
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options] <dataset>"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python (optional).")
    parser.add_option(
        "", "--bfrm_bin", dest="bfrm_bin", default=None,
        help="Specify the path to the BFRM binary.")
    parser.add_option(
        "", "--arrayplot", dest="arrayplot", default="arrayplot.py",
        help="Specify the command to run arrayplot.")
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
        "-z", "--archive", dest="archive", action="store_true", default=None,
        help="Archive the raw output.  Helpful for GenePattern.")
    
    group = OptionGroup(parser, "Filtering")
    group.add_option(
        "--filter_mean", dest="filter_mean", type=float, default=None,
        help="Remove this portion of genes based on mean expression.")
    group.add_option(
        "--filter_var", dest="filter_var", type=float, default=None,
        help="Remove this portion of genes based on variance.")
    group.add_option(
        "--cutoff", dest="cutoff", type=float, default=0.99,
        help="Cutoff probability for a gene to be in a factor.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "BFRM Parameters")
    group.add_option(
        "--nc", dest="num_control_vars", type="int", default=None,
        help="Specify the number of control variables to use.")
    group.add_option(
        "--start_factors", dest="start_factors", type="int", default=None,
        help="The number of factors to start the analysis.")
    group.add_option(
        "--nucleus_file", dest="nucleus_file", default=None,
        help="A file that contains the genes to start the evolution.")
    group.add_option(
        "--evol_max_factors", dest="evol_max_factors", default=None,
        help="Maximum number of factors for the evolution.")
    group.add_option(
        "--evol_max_genes", dest="evol_max_genes", default=None,
        help="Maximum number of genes for the evolution.")
    parser.add_option_group(group)

    # Parse the arguments.
    options, args = parser.parse_args()

    if options.cutoff <= 0 or options.cutoff > 1:
        parser.error("Cutoff probability should be between 0 and 1.")
    if options.filter_mean and (
        options.filter_mean < 0 or options.filter_mean >= 1):
        parser.error("filter_mean filter should be between 0 and 1.")
    if options.filter_var and (
        options.filter_var < 0 or options.filter_var >= 1):
        parser.error("filter_var filter should be between 0 and 1.")

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import this after the library path is set.
    import arrayio
    from genomicode import archive
    from genomicode import genepattern

    genepattern.fix_environ_path()

    if len(args) != 1:
        parser.error("Please specify a file to factor.")
    filename, = args
    assert os.path.exists(filename), "File not found: %s" % filename

    if options.nucleus_file:
        assert os.path.exists(options.nucleus_file), "I could not find: %s" % \
               options.nucleus_file

    # Set up the files.
    file_layout = make_file_layout(options.outpath)
    init_paths(file_layout)

    # Read the matrix and convert to GCT format.
    x = arrayio.read(filename)
    MATRIX_orig = arrayio.convert(x, to_format=arrayio.gct_format)
    print "Read data set with %d genes and %d samples." % (
        MATRIX_orig.nrow(), MATRIX_orig.ncol())

    # Make a copy so that in-place changes (like log_matrix) won't
    # affect the original matrix.
    MATRIX = MATRIX_orig.matrix()

    # Log the data set if necessary.
    log_matrix(MATRIX)

    # Filter out based on mean and varian
    MATRIX = filter_dataset(MATRIX, options.filter_mean, options.filter_var)
    if MATRIX.nrow() != MATRIX_orig.nrow():
        print "Filtered from %d genes to %d." % (
            MATRIX_orig.nrow(), MATRIX.nrow())
    
    # Write out the data sets.
    write_dataset(file_layout.DATASET_ORIG, MATRIX_orig)
    write_dataset(file_layout.DATASET, MATRIX)

    # Run BFRM.
    run_bfrm(
        file_layout, options.bfrm_bin, options.num_control_vars,
        options.start_factors, options.nucleus_file, options.evol_max_factors,
        options.evol_max_genes)

    # Generate output files.
    summarize_factor_scores(
        file_layout, options.cutoff, options.python, options.arrayplot,
        options.cluster, options.libpath)
    summarize_gene_factor_probs(
        file_layout, options.cutoff, options.python, options.arrayplot,
        options.cluster, options.libpath)
    summarize_factor_geneset(file_layout, options.cutoff)

    # BFRM model file should always be archived.
    archive.zip_path(file_layout.BFRM, noclobber=False)
    
    if options.archive:
        print "Archiving results."
        archive.zip_path(file_layout.ATTIC, noclobber=False)

    print "Done."

    
if __name__ == '__main__':
    main()
