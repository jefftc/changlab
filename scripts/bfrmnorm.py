#!/usr/bin/env python

import os, sys

# make_file_layout
# init_paths
# 
# read_matrices
# log_matrices
# write_dataset
#
# run_bfrm
# 
# summarize_dataset
# summarize_filtered_genes
# summarize_heatmaps
# summarize_pca
# summarize_report

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
            File.DS_ORIG("dataset.original.gct"),
            File.DS_PROC_FILTERED("pre_norm.filtered.gct"),
            File.DS_PROC_COORD("pre_norm.pca.txt"),
            File.DS_PROC_POV("pre_norm.pov"),
            File.DS_FINAL_FILTERED("normalized.filtered.gct"),
            File.DS_FINAL_COORD("normalized.pca.txt"),
            File.DS_FINAL_POV("normalized.pov"),
            ),
        Path.BFRM("bfrm",
            File.BFRM_PROBE_IDS("probeids.txt"),
            File.BFRM_SAMPLE_IDS("sids.txt"),
            File.BFRM_DATASET("dataset.txt"),
            File.BFRM_PARAMETERS("parameters.txt"),
            File.BFRM_CORRECTED("correctedData.txt"),
            File.BFRM_MA("mA.txt"),
            ),
        
        File.DS_PROC("pre_norm.gct"),
        File.DS_PROC_HEATMAP("pre_norm.heatmap.png"),
        File.DS_PROC_SCATTER("pre_norm.pca.png"),
        
        File.DS_FINAL("normalized.gct"),
        File.DS_FINAL_HEATMAP("normalized.heatmap.png"),
        File.DS_FINAL_SCATTER("normalized.pca.png"),
        
        File.REPORT("REPORT.html"),
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

def read_matrices(filenames):
    from genomicode import binreg

    if not filenames:
        return []

    x = binreg.read_matrices(filenames)
    DATA, ALIGNED = x

    for d, filename in zip(DATA, filenames):
        f = os.path.split(filename)[1]
        print "%s has %d genes and %d samples." % (f, d.nrow(), d.ncol())
    print "Merged file has %d genes." % ALIGNED[0].nrow()
    sys.stdout.flush()

    return ALIGNED

def log_matrices(names, matrices):
    # Log each variable if necessary.  Will log in place.  Return a
    # boolean indicating whether anything was logged.  test can be None.
    from genomicode import jmath
    from genomicode import binreg

    any_files_logged = False
    for name, matrix in zip(names, matrices):
        msg = "I will not log %s." % name

        if not binreg.is_logged_array_data(matrix):
            msg = "I will log %s." % name
            matrix._X = jmath.log(matrix._X, base=2, safe=1)
            any_files_logged = True
        print msg
    sys.stdout.flush()
    return any_files_logged

def write_dataset(filename, matrices):
    import arrayio
    from genomicode import binreg

    DATA = binreg.merge_gct_matrices(*matrices)
    arrayio.gct_format.write(DATA, open(filename, 'w'))

def run_bfrm(bfrm_path, num_factors, file_layout, matlab_bin):
    import re
    import subprocess
    import arrayio
    from genomicode import filelib
    from genomicode import bfrm
    from genomicode import matlab

    assert filelib.exists_nz(file_layout.DS_PROC)
    DATA = arrayio.read(file_layout.DS_PROC)

    # Write the data to a file.
    data = DATA.value()
    handle = open(file_layout.BFRM_DATASET, 'w')
    for x in data:
        print >>handle, "\t".join(map(str, x))
    handle.close()

    # Write the probeids.txt and sampleids.txt files.
    sample_ids = DATA.col_names(arrayio.COL_ID)
    probe_ids = DATA.row_names("NAME")
    x = ["%s\n" % x for x in sample_ids]
    open(file_layout.BFRM_SAMPLE_IDS, 'w').writelines(x)
    x = ["%s\n" % x for x in probe_ids]
    open(file_layout.BFRM_PROBE_IDS, 'w').writelines(x)

    bfrm_path = bfrm.find_bfrm_normalize(bfrm_path)
    assert bfrm_path is not None, "Could not find BFRM normalize code."
    assert os.path.exists(bfrm_path)

    # Use Matlab to set up the files.
    print "Initializing files for BFRM."
    PARAMETERS = [
        ("root", "'%s'" % bfrm_path),
        ("NUM_CONTROL_FACTORS", num_factors),
        ]
    setup_file = os.path.join(bfrm_path, "setup.m")
    assert filelib.exists_nz(setup_file)
    lines = open(setup_file).readlines()
    for key, value in PARAMETERS:
        # Find the key in the parameters file.
        for i, line in enumerate(lines):
            if re.match(r"%s\s*=" % key, line, re.IGNORECASE):
                break
        else:
            raise AssertionError, "I could not find parameter: %s" % key
        lines[i] = "%s = %s;\n" % (key, value)

    # Run setup.m from Matlab.
    script = "".join(lines)
    x = matlab.run(
        script, matlab_bin=matlab_bin, working_path=file_layout.BFRM)
    print x
    assert filelib.exists_nz(file_layout.BFRM_PARAMETERS)

    # Execute BFRM.
    print "BFRM normalizing with %d factors." % num_factors
    sys.stdout.flush()
    bfrm_bin = os.path.join(bfrm_path, "bfrm64")
    assert filelib.exists_nz(bfrm_bin)
    cwd = os.getcwd()
    try:
        os.chdir(file_layout.BFRM)
        cmd = "%s %s" % (bfrm_bin, file_layout.BFRM_PARAMETERS)
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
    assert filelib.exists_nz(file_layout.BFRM_MA)

    # Run computeCorrected.m from Matlab.
    correct_file = os.path.join(bfrm_path, "computeCorrected.m")
    assert filelib.exists_nz(correct_file)
    script = open(correct_file).read()
    x = matlab.run(
        script, matlab_bin=matlab_bin, working_path=file_layout.BFRM)
    print x
    assert filelib.exists_nz(file_layout.BFRM_CORRECTED)

def summarize_dataset(file_layout):
    import arrayio
    from genomicode import Matrix
    from genomicode import filelib

    DATA_orig = arrayio.read(file_layout.DS_ORIG)

    x = arrayio.tab_delimited_format.read(file_layout.BFRM_CORRECTED)
    M_corrected = x._X
    
    assert M_corrected
    nrow, ncol = len(M_corrected), len(M_corrected[0])
    assert (nrow, ncol) == DATA_orig.dim(), "%s %s" % (
        (nrow, ncol), DATA_orig.dim())

    DATA_corrected = DATA_orig.matrix()
    DATA_corrected._X = M_corrected
    handle = open(file_layout.DS_FINAL, 'w')
    arrayio.gct_format.write(DATA_corrected, handle)
    handle.close()

def summarize_filtered_genes(file_layout):
    # Select the <NUM_GENES> genes that vary most by variance.
    import arrayio
    from genomicode import binreg
    from genomicode import pcalib

    # If this is changed, you also need to change the report.
    NUM_GENES = 1000

    DATA_orig = arrayio.read(file_layout.DS_PROC)
    DATA_final = arrayio.read(file_layout.DS_FINAL)
    assert binreg.are_rows_aligned(DATA_orig, DATA_final)

    # Select the genes with the greatest variance.
    I = pcalib.select_genes_var(DATA_orig._X, NUM_GENES)
    DATA_orig = DATA_orig.matrix(I, None)
    DATA_final = DATA_final.matrix(I, None)
    
    arrayio.gct_format.write(
        DATA_orig, open(file_layout.DS_PROC_FILTERED, 'w'))
    arrayio.gct_format.write(
        DATA_final, open(file_layout.DS_FINAL_FILTERED, 'w'))

def summarize_heatmaps(
    python, arrayplot, cluster, file_layout, libpath=[]):
    import arrayio
    from genomicode import plotlib

    # Load the data sets.
    DATA_orig = arrayio.gct_format.read(file_layout.DS_PROC_FILTERED)
    DATA_final = arrayio.gct_format.read(file_layout.DS_FINAL_FILTERED)
    assert DATA_final.dim() == DATA_orig.dim()

    nrow, ncol = DATA_orig.dim()
    x = plotlib.find_tall_heatmap_size(
        #nrow, ncol, min_box_width=20, max_total_height=2000,
        nrow, ncol, max_total_height=2000, max_total_width=2000)
    xpix, ypix = x
    #print "SIZE", nrow, ncol, xpix, ypix
    
    plotlib.plot_heatmap(
        file_layout.DS_FINAL_FILTERED, file_layout.DS_FINAL_HEATMAP,
        xpix, ypix, color="bild", show_colorbar=True, show_grid=True,
        gene_center="mean", gene_normalize="var",
        array_label=True, 
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)
    plotlib.plot_heatmap(
        file_layout.DS_PROC_FILTERED, file_layout.DS_PROC_HEATMAP,
        xpix, ypix, color="bild", show_colorbar=True, show_grid=True,
        gene_center="mean", gene_normalize="var",
        array_label=True, 
        python=python, arrayplot=arrayplot, cluster=cluster, libpath=libpath)

def _make_scatter(povray, pca_file, pov_file, out_file):
    from genomicode import filelib
    from genomicode import pcalib

    ds = [d for d in filelib.read_row(pca_file, header=1)]
    X = [float(d.PC_0) for d in ds]
    Y = [float(d.PC_1) for d in ds]
    Z = [float(d.PC_2) for d in ds]
    DATASET = [int(d.Dataset) for d in ds]
    assert min(DATASET) >= 0 and max(DATASET) < 256

    x = pcalib.plot_scatter(
        X, Y, out_file, group=DATASET, pov_file=pov_file, povray=povray)
    print x
    sys.stdout.flush()
    assert os.path.exists(out_file), "Failed to plot predictions."

def summarize_pca(povray, file_layout, matrices):
    import arrayio
    from genomicode import pcalib

    # Load the data sets.
    DATA_orig = arrayio.gct_format.read(file_layout.DS_PROC_FILTERED)
    DATA_final = arrayio.gct_format.read(file_layout.DS_FINAL_FILTERED)
    assert DATA_final.dim() == DATA_orig.dim()
    nrow, ncol = DATA_orig.dim()
    nmin = min(nrow, ncol)

    # Get the names of the samples.
    samples = []
    for m in matrices:
        x = m.col_names(arrayio.COL_ID)
        samples.extend(x)
    assert samples == DATA_orig.col_names(arrayio.COL_ID)
    assert samples == DATA_final.col_names(arrayio.COL_ID)

    # Get assignment of data set to sample.
    dataset = []
    for i, m in enumerate(matrices):
        x = [i]*m.ncol()
        dataset.extend(x)
    assert len(dataset) == len(samples)

    FL = file_layout
    K = 3
    PC_orig = pcalib.svd_project_cols(DATA_orig._X, K)
    PC_final = pcalib.svd_project_cols(DATA_final._X, K)
    _write_svd_coord(FL.DS_PROC_COORD, PC_orig, samples, dataset)
    _write_svd_coord(FL.DS_FINAL_COORD, PC_final, samples, dataset)
    _make_scatter(povray, FL.DS_PROC_COORD, FL.DS_PROC_POV, FL.DS_PROC_SCATTER)
    _make_scatter(
        povray, FL.DS_FINAL_COORD, FL.DS_FINAL_POV, FL.DS_FINAL_SCATTER)

def summarize_report(
    filenames, matrices, num_factors, start_time, file_layout):
    import time
    import subprocess
    from genomicode import htmllib
    from genomicode import parselib

    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")

    assert len(filenames) == len(matrices)

    lines = []
    w = lines.append
    w("<HTML>")
    w(htmllib.HEAD(htmllib.TITLE("BFRMNormalize Report")))
    w("<BODY>")
    w(htmllib.CENTER(htmllib.H1(
        htmllib.EM("BFRMNormalize") + " Report")))

    w(htmllib.H3("I.  Overview"))
    l = []
    l.append("I normalized the following data sets with %d factors." %
             num_factors)
    files = [os.path.split(x)[1] for x in filenames]
    for i in range(len(files)):
        name = files[i]
        num_samples = matrices[i].ncol()
        x = "%s (%d samples)" % (name, num_samples)
        l.append(htmllib.LI() + x)
    l = "\n".join(l)
    w(htmllib.UL(l))
    
    w(htmllib.P())
    x = os.path.split(file_layout.DS_PROC)[1]
    w("The merged gene expression data set is available at " +
      htmllib.A(x, href=x) + ".")
    w(htmllib.BR())
    x = os.path.split(file_layout.DS_FINAL)[1]
    w("The normalized data set is available at " +
      htmllib.A(x, href=x) + ".")
    
    w(htmllib.P())
    w(htmllib.H3("II.  Results"))

    # Make the table of the heatmaps.
    x = os.path.split(file_layout.DS_PROC_HEATMAP)[1]
    x1 = htmllib.CENTER(
        htmllib.B("Before Normalization") + htmllib.BR() +
        htmllib.A(htmllib.IMG(height=480, src=x), href=x)
        )
    x = os.path.split(file_layout.DS_FINAL_HEATMAP)[1]
    x2 = htmllib.CENTER(
        htmllib.B("After Normalization") + htmllib.BR() +
        htmllib.A(htmllib.IMG(height=480, src=x), href=x)
        )
    row1 = htmllib.TR(htmllib.TD(x1) + htmllib.TD(x2))

    # Make sure this is really the 1000 genes.
    x = htmllib.TD(
        htmllib.B("Figure 1: Heatmaps. ") +
        "These heatmaps show the expression patterns in the data before "
        "and after normalization.  "
        "The rows contain the 1000 genes that exhibit the highest variance "
        "in gene expression profile across the original data set.  "
        "The columns contain the samples in the data sets provided.  "
        "The genes and samples are in the same order in both heatmaps.  "
        "Warm colors indicate high expression of the gene, and cool colors "
        "indicate low expression.",
        colspan=2)
    row2 = htmllib.TR(x)

    w(htmllib.TABLE(row1 + row2, border=0, cellspacing=10, width="50%%"))
    
    w(htmllib.P())

    # Make the table of the scatter plots.
    x = os.path.split(file_layout.DS_PROC_SCATTER)[1]
    x1 = htmllib.CENTER(
        htmllib.B("Before Normalization") + htmllib.BR() +
        htmllib.A(htmllib.IMG(height=400, src=x), href=x)
        )
    x = os.path.split(file_layout.DS_FINAL_SCATTER)[1]
    x2 = htmllib.CENTER(
        htmllib.B("After Normalization") + htmllib.BR() +
        htmllib.A(htmllib.IMG(height=400, src=x), href=x)
        )
    row1 = htmllib.TR(htmllib.TD(x1) + htmllib.TD(x2))

    x = htmllib.TD(
        htmllib.B("Figure 2: PCA Plots. ") +
        "These plots show the gene expression profiles of the samples "
        "plotted on the first two principal components.  "
        "Each point represents a sample, and samples from the same data "
        "set have the same color.  "
        "If there are batch effects, the samples from the same data set "
        "(the same color) will cluster together.  "
        "If there are no batch effects, the colors should be mixed.",
        colspan=2)
    row2 = htmllib.TR(x)

    w(htmllib.TABLE(row1 + row2, border=0, cellspacing=10, width="50%%"))
    

    # Format the current time.
    end_time = time.time()
    time_str = parselib.pretty_date(start_time)
    x = int(end_time-start_time)
    num_min = x / 60
    num_secs = x % 60
    if num_min == 0:
        run_time = "%ss" % parselib.pretty_int(num_secs)
    else:
        run_time = "%sm %ss" % (parselib.pretty_int(num_min), num_secs)

    # Get the hostname.
    cmd = "hostname"
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    wh, r = p.stdin, p.stdout
    wh.close()
    hostname = r.read().strip()
    assert hostname, "I could not get the hostname."

    w(htmllib.P())
    w(htmllib.HR())
    w(htmllib.EM(
        "This analysis was run on %s on %s.  It took %s to complete." %
        (time_str, hostname, run_time)))

    w("</BODY>")
    w("</HTML>")

    x = "\n".join(lines) + "\n"
    outfile = file_layout.REPORT
    open(outfile, 'w').write(x)

def _write_svd_coord(filename, PC, samples, dataset):
    if not PC:
        return
    
    x = ["PC-%d" % x for x in range(len(PC[0]))]
    header = ["SampleID", "Dataset"] + x
    
    handle = open(filename, 'w')
    print >>handle, "\t".join(header)
    for s, d, xyz in zip(samples, dataset, PC):
        x = [s, d] + xyz
        print >>handle, "\t".join(map(str, x))
    handle.close()
    
def main():
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options] <file1> <file2> ..."
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "-f", "--num_factors", dest="num_factors", type="int", 
        default=15, help="Number of factors to use for normalization.")
    parser.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python (optional).")
    parser.add_option(
        "", "--bfrm", dest="bfrm_path", default=None,
        help="Specify the path to the BFRM_normalize directory.")
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
        "-z", "--archive", dest="archive", action="store_true", default=None,
        help="Archive the raw output.  Helpful for GenePattern.")
    

    # Parse the arguments.
    options, args = parser.parse_args()

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import this after the library path is set.
    import time
    import arrayio
    from genomicode import filelib
    from genomicode import archive
    from genomicode import genepattern

    start_time = time.time()

    genepattern.fix_environ_path()

    if not args:
        parser.error("Please specify files to normalize.")
    filenames = args
    names = [os.path.split(x)[-1] for x in filenames]
    for filename in filenames:
        assert filelib.exists(filename), "File not found: %s" % filename

    # Check to make sure value for num_factors is reasonable.
    MIN_FACTORS, MAX_FACTORS = 2, 100
    if options.num_factors < MIN_FACTORS:
        parser.error("At least %d factors are required." % MIN_FACTORS)
    elif options.num_factors > MAX_FACTORS:
        parser.error("%d factors is too many.  Maximum is %d." % (
            options.num_factors, MAX_FACTORS))

    # Set up the files.
    file_layout = make_file_layout(options.outpath)
    init_paths(file_layout)

    # Read each of the input files and align them.
    matrices = read_matrices(filenames)

    # Make sure the number of factors don't exceed the size of the
    # matrices.
    if matrices and options.num_factors > matrices[0].nrow():
        parser.error("Too many factors.")

    # Standardize each of the matrices to GCT format.
    for i in range(len(matrices)):
        matrices[i] = arrayio.convert(
            matrices[i], to_format=arrayio.gct_format)
    write_dataset(file_layout.DS_ORIG, matrices)

    # Log each of the matrices if needed.
    log_matrices(names, matrices)
    write_dataset(file_layout.DS_PROC, matrices)
    sys.stdout.flush()

    # Format the parameters and output files for bfrm.
    run_bfrm(
        options.bfrm_path, options.num_factors, file_layout, options.matlab)

    # Generate some files for output.
    summarize_dataset(file_layout)
    summarize_filtered_genes(file_layout)
    summarize_heatmaps(
        options.python, options.arrayplot, options.cluster,
        file_layout, options.libpath)
    summarize_pca(options.povray, file_layout, matrices)
    summarize_report(
        filenames, matrices, options.num_factors, start_time, file_layout)

    # Archive the BFRM stuff, and the big files.
    if options.archive:
        print "Archiving results."
        archive.zip_path(file_layout.BFRM, noclobber=False)
        archive.zip_path(file_layout.ATTIC, noclobber=False)
        #archive.zip_path(file_layout.DS_PROC, noclobber=False)
        #archive.zip_path(file_layout.DS_FINAL, noclobber=False)

    print "Done."

    
if __name__ == '__main__':
    main()
