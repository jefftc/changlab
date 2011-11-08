#!/usr/bin/env python

import os, sys

# make_file_layout
# init_paths
#
# write_dataset
# write_model
# write_bfrm_dataset
# write_sample_probe_ids
# write_bfrm_files
# check_model
#
# run_bfrm_project
# log_matrix
#
# summarize_factor_scores

def make_file_layout(outpath):
    from genomicode import filelayout

    outpath = outpath or "."
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath = os.path.realpath(outpath)

    Path, File = filelayout.Path, filelayout.File
    
    # Have only set set of these files for the whole analysis.
    FILES = [
        Path.ATTIC("attic"),
        File.DATASET("dataset.gct"),
        File.FACTOR_SCORES("factors.pcl"),
        File.FACTOR_CDT("factors.cdt"),   # created by cluster
        File.FACTOR_ATR("factors.atr"),   # created by cluster
        File.FACTOR_GTR("factors.gtr"),   # created by cluster
        File.FACTOR_SCORES_PNG("factors.png"),
        Path.BFRM("bfrm",
            File.BFRM_DATASET("dataset.txt"),
            File.BFRM_SPROBE_IDS("probeidsSmp.txt"),
            File.BFRM_MA("mA.txt"),
            File.BFRM_MPOSTPIB("mPostPib.txt"),
            File.BFRM_MPSI("mPsi.txt"),
            File.BFRM_MVARIABLESIN("mVariablesIn.txt"),  # optional
            File.BFRM_PROBE_IDS("probeids.txt"),
            File.BFRM_AF("af.txt"),  # sample x factor
            File.BFRM_Y("Y.txt"),    # gene x sample
            ),
        File.BFRM_MODEL("bfrm_model.zip"),
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

def write_model(filename, file_layout):
    check_model(filename)
    x = open(filename).read()
    open(file_layout.BFRM_MODEL, 'w').write(x)

def write_bfrm_dataset(filename, DATA):
    data = DATA.value()
    handle = open(filename, 'w')
    for x in data:
        print >>handle, "\t".join(map(str, x))
    handle.close()

def write_sample_probe_ids(filename, DATA):
    from genomicode import bfrm
    
    name = bfrm.get_affy_row_name(DATA)
    probeset_ids = DATA.row_names(name)
    x = ["%s\n" % x for x in probeset_ids]
    open(filename, 'w').writelines(x)

def write_bfrm_files(path, model_file):
    import zipfile
    from genomicode import archive

    opj = os.path.join
    
    zfile = zipfile.ZipFile(model_file)
    s2f = archive.unzip_dict(model_file)

    # list of filename, required
    files = [
        ("mA.txt", 1),
        ("mPostPib.txt", 1),
        ("mPsi.txt", 1),
        ("mVariablesIn.txt", 0),
        ("probeids.txt", 1),
        ]
    for name, required in files:
        if name not in s2f and not required:
            continue
        assert name in s2f, "I could not find '%s' in the model." % name
        x = zfile.open(s2f[name]).read()
        open(opj(path, name), 'w').write(x)

def check_model(filename):
    from genomicode import archive
    
    s2f = archive.unzip_dict(filename)
    assert "mA.txt" in s2f
    assert "mPostPib.txt" in s2f
    assert "mPsi.txt" in s2f
    assert "probeids.txt" in s2f

def run_bfrm_project(file_layout, bfrm_path, matlab_bin):
    import arrayio
    from genomicode import bfrm
    from genomicode import matlab

    param_file = "parameters.txt"
    model = bfrm.read_clean_model(
        file_layout.BFRM_MODEL, param_file=param_file)
    num_factors = len(model["FACTOR_O"])
    assert num_factors, "No latent factors in the BFRM model."
    x = "Projecting %d latent factors onto data set." % num_factors
    if num_factors == 1:
        x = x.replace("factors", "factor")
    print x
    
    DATA = arrayio.read(file_layout.DATASET)
    
    bfrm_path = bfrm.find_bfrm_project(bfrm_path)
    assert bfrm_path is not None, "I could not find BFRM_project."
    bfrm_path = os.path.realpath(bfrm_path)

    # Write out the dataset and probe IDs.
    write_bfrm_dataset(file_layout.BFRM_DATASET, DATA)
    write_sample_probe_ids(file_layout.BFRM_SPROBE_IDS, DATA)

    # Write the BFRM model files.
    write_bfrm_files(file_layout.BFRM, file_layout.BFRM_MODEL)

    # Make sure some of the probes are the same.
    pid = [x.strip() for x in open(file_layout.BFRM_PROBE_IDS)]
    pid = [pid[i] for i in model["VariablesIn"]]
    spid = [x.strip() for x in open(file_layout.BFRM_SPROBE_IDS)]
    pid = [x.lower() for x in pid]
    spid = [x.lower() for x in spid]
    intersect = [x for x in pid if x in spid]
    assert intersect, "No common probes between model and data set."
    if len(intersect) < len(pid):
        x = "Warning: model contains %d probe IDs, but only matched " + \
            "%d in data set."
        print x % (len(pid), len(intersect))

    # Run the matlab script.
    lines = []
    w = lines.append
    w("addpath '%s';\n" % bfrm_path)
    w("addpath '%s/bfrm';\n" % bfrm_path)
    w("y = load('%s');\n" % file_layout.BFRM_DATASET)
    w("probeidsSmp = readWordlist('%s');\n" % file_layout.BFRM_SPROBE_IDS)
    w("[af Y sampleids] = getFacScores('%s/', y, probeidsSmp);" %
      file_layout.BFRM)
    w("save('%s', 'af', '-ASCII', '-TABS');\n" % file_layout.BFRM_AF)
    w("save('%s', 'Y', '-ASCII', '-TABS');\n" % file_layout.BFRM_Y)
    script = "".join(lines)
    x = matlab.run(
        script, matlab_bin=matlab_bin, working_path=file_layout.OUTPATH)
    print x
    sys.stdout.flush()

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

def summarize_factor_scores(file_layout, python, arrayplot, cluster, libpath):
    import zipfile
    import arrayio
    from genomicode import Matrix
    from genomicode import jmath
    from genomicode import archive
    from genomicode import graphlib
    from genomicode import bfrm

    DATA = arrayio.read(file_layout.DATASET)

    param_file = "parameters.txt"
    model = bfrm.read_clean_model(
        file_layout.BFRM_MODEL, param_file=param_file)
    num_factors = model["F"].nrow()

    # Load the factor names.
    assert zipfile.is_zipfile(file_layout.BFRM_MODEL)
    s2f = archive.unzip_dict(file_layout.BFRM_MODEL)
    assert "factorids.txt" in s2f, "Missing: factorids.txt"
    zfile = zipfile.ZipFile(file_layout.BFRM_MODEL)
    factor_names = [x.strip() for x in zfile.open(s2f["factorids.txt"])]
    assert len(factor_names) == num_factors
    
    # sample x factor matrix
    F = arrayio.read(file_layout.BFRM_AF)
    assert F.nrow() == DATA.ncol()
    F_X = jmath.transpose(F._X)

    # F_X contains all factors, including intercept and design.
    # Remove all but the latent factors.
    F_X = F_X[-num_factors:]

    # Sort the factors so they'll be in the same order as the clean
    # model.
    assert len(F_X) == len(model["FACTOR_O"])
    F_X = [F_X[i] for i in model["FACTOR_O"]]
    factor_names = [factor_names[i] for i in model["FACTOR_O"]]

    # Write out the projected factor scores.
    SAMPLE_NAME = arrayio.tdf.SAMPLE_NAME
    row_names = {}
    col_names = {}
    row_names["xID"] = factor_names
    col_names[SAMPLE_NAME] = DATA.col_names(SAMPLE_NAME)
    M = Matrix.InMemoryMatrix(F_X, row_names, col_names)
    arrayio.pcl_format.write(M, file_layout.FACTOR_SCORES)

    # Make the heatmap.
    x = graphlib.find_wide_heatmap_size(
        M.nrow(), M.ncol(), min_box_height=10, min_box_width=10,
        max_total_height=768, max_total_width=1024)
    xpix, ypix = x
    ypix = min(ypix, xpix*4)
    x = graphlib.plot_heatmap(
        file_layout.FACTOR_SCORES, file_layout.FACTOR_SCORES_PNG,
        xpix, ypix,
        color="bild", show_colorbar=True, show_grid=True,
        gene_center="mean", gene_normalize="var",
        gene_label=True, cluster_genes=True,
        array_label=True, cluster_arrays=True,
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

    # Clean up the cluster files.
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

def main():
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options] <bfrm_model> <dataset>"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "", "--bfrm_path", dest="bfrm_path", default=None,
        help="Specify the path to BFRM_project.")
    parser.add_option(
        "", "--matlab", dest="matlab", default="matlab",
        help="Specify the command to run matlab.")
    parser.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python (optional).")
    parser.add_option(
        "", "--arrayplot", dest="arrayplot", default=None,
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
    
    # Parse the arguments.
    options, args = parser.parse_args()

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import this after the library path is set.
    import arrayio
    from genomicode import archive
    from genomicode import genepattern

    genepattern.fix_environ_path()

    if len(args) != 2:
        parser.error("Please specify files.")
    model_file, filename = args
    assert os.path.exists(model_file), "File not found: %s" % model_file
    assert os.path.exists(filename), "File not found: %s" % filename

    # Set up the files.
    file_layout = make_file_layout(options.outpath)
    init_paths(file_layout)

    # Read the matrix and convert to GCT format.
    x = arrayio.read(filename)
    MATRIX = arrayio.convert(x, to_format=arrayio.gct_format)
    print "Read data set with %d genes and %d samples." % (
        MATRIX.nrow(), MATRIX.ncol())

    log_matrix(MATRIX)

    # Write out the data sets.
    write_dataset(file_layout.DATASET, MATRIX)

    # Save the BFRM model.
    write_model(model_file, file_layout)

    # Run BFRM projection.
    run_bfrm_project(
        file_layout, options.bfrm_path, options.matlab)

    # Generate output files.
    summarize_factor_scores(
        file_layout, options.python, options.arrayplot, options.cluster,
        options.libpath)

    if options.archive:
        print "Archiving results."
        archive.zip_path(file_layout.ATTIC, noclobber=False)
        archive.zip_path(file_layout.BFRM, noclobber=False)

    print "Done."

    
if __name__ == '__main__':
    main()
