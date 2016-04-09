#!/usr/bin/env python


def choose_gene_names(MATRIX):
    # Return tuple of (header for gene_id, header for gene_names).
    # Either of the headers can be None.
    from genomicode import arrayplatformlib as apl

    if not MATRIX.row_names():
        return None, None
    if len(MATRIX.row_names()) == 1:
        return MATRIX.row_names()[0], None
    if len(MATRIX.row_names()) == 2:
        return MATRIX.row_names()

    # By default, set as the first two rows.
    geneid_header, genename_header = MATRIX.row_names()[:2]

    # If I can find a better header for the annotation, then use it.
    x = apl.score_matrix(MATRIX, min_score=0.5)
    for score in x:
        # TODO: Use categories.
        if score.platform_name in [
            "Entrez_Symbol_human", "Entrez_Symbol_mouse"]:
            genename_header = score.header
    return geneid_header, genename_header


def _write_matrix(outhandle, header, DATA_py):
    if type(outhandle) is type(""):
        outhandle = open(outhandle, 'w')

    print >>outhandle, "\t".join(header)
    for x in DATA_py:
        assert len(x) == len(header)
        # Convert None to "".
        for i in range(len(x)):
            if x[i] is None:
                x[i] = ""
        print >>outhandle, "\t".join(map(str, x))
    outhandle.flush()
    
    
def find_diffexp_genes(
    outfile, gmt_file, algorithm, paired,
    MATRIX, geneid_header, genename_header, genename_delim,
    name1, name2, classes, fold_change, p_cutoff, fdr_cutoff, bonf_cutoff, 
    sam_DELTA, sam_qq_file, edger_tagwise_dispersion, num_procs):
    # classes must be 0, 1, None.
    import os
    import sys
    import math
    import StringIO
    import warnings
    
    from rpy2 import rinterface
    
    from genomicode import config
    from genomicode import jmath
    from genomicode import genesetlib

    algorithm2function_unpaired = {
        "fold_change" : "find.de.genes.fc",
        "ttest" : "find.de.genes.ttest",
        "sam" : "find.de.genes.sam",
        "ebayes" : "find.de.genes.ebayes",
        "deseq2" : "find.de.genes.deseq2",
        "edger" : "find.de.genes.edgeR",
        }
    algorithm2function_paired = {
        "ebayes" : "find.de.genes.paired.ebayes",
        }
    algorithm2function = algorithm2function_unpaired
    if paired:
        algorithm2function = algorithm2function_paired
        assert algorithm in algorithm2function_paired, \
               "No paired version of %s" % algorithm
    assert algorithm in algorithm2function, "Unknown algorithm: %s" % algorithm

    
    # Select the relevant columns from MATRIX.
    I = [i for (i, x) in enumerate(classes) if x in [0, 1]]
    assert len(I)
    MATRIX = MATRIX.matrix(None, I)
    classes = [classes[i] for i in I]

    # All algorithms except "fold_change" need at least 2 samples of
    # each class.
    counts = {}
    for x in classes:
        counts[x] = counts.get(x, 0) + 1
    assert sorted(counts) == [0, 1], "Only one class represented."
    if algorithm != "fold_change":
        assert counts[0] >= 2, "There must be at least 2 of each class."
        assert counts[1] >= 2, "There must be at least 2 of each class."

    names = [name1, name2]
    X = MATRIX._X
    Y = [names[x] for x in classes]
    sample_name = None
    if MATRIX.col_names():
        sample_name = MATRIX.col_names(MATRIX.col_names()[0])

    x = choose_gene_names(MATRIX)
    if not geneid_header:
        geneid_header = x[0]
    if not genename_header:
        genename_header = x[1]
    assert not geneid_header or geneid_header in MATRIX.row_names()
    assert not genename_header or genename_header in MATRIX.row_names()

    R = jmath.start_R()
    de_lib = os.path.join(config.changlab_Rlib, "diffexp.R")
    stat_lib = os.path.join(config.changlab_Rlib, "statlib.R")
    assert os.path.exists(de_lib), "I could not find file: %s" % de_lib
    assert os.path.exists(stat_lib), "I could not find file: %s" % stat_lib
    R('source("%s")' % de_lib)
    R('source("%s")' % stat_lib)

    jmath.R_equals(X, "X")
    jmath.R_equals(Y, "Y")
    if sample_name:
        jmath.R_equals(sample_name, "sample.name")
        jmath.R('colnames(X) <- sample.name')

    geneid = genenames = None
    if geneid_header:
        geneid = MATRIX.row_names(geneid_header)
        jmath.R_equals(geneid, "geneid")
    if genename_header:
        genenames = MATRIX.row_names(genename_header)
        jmath.R_equals(genenames, "genenames")

    # Set up the arguments.
    args = ["X", "Y"]
    if algorithm == "sam":
        args.append("%g" % sam_DELTA)
    if geneid:
        args.append("geneid=geneid")
    if genenames:
        args.append("genenames=genenames")
    # Pass the fold change to the algorithm, because it can affect the
    # multiple hypothesis correction.
    if fold_change is not None:
        args.append("FOLD.CHANGE=%g" % fold_change)
    if algorithm in ["ttest", "deseq2"]:
        args.append("NPROCS=%d" % num_procs)  # t-test only
    #if show_all_genes and algorithm != "sam":
    if algorithm not in ["sam", "fold_change"]:
        args.append("filter.p05=FALSE")
    if algorithm == "edger":
        if edger_tagwise_dispersion:
            args.append("tagwise.dispersion=TRUE")
        else:
            args.append("tagwise.dispersion=FALSE")

    # Prevent SAM from writing junk to the screen.
    handle = StringIO.StringIO()
    old_stdout = sys.stdout
    sys.stdout = handle

    # Call the proper R function.  DESeq2 throws off a lot of
    # warnings.  Turn them off temporarily.
    fn = algorithm2function[algorithm]
    x = ", ".join(args)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        R("x <- %s(%s)" % (fn, x))
    R("DATA <- x$DATA")
    DATA_R = R["DATA"]

    sys.stdout = old_stdout

    # Write out a QQ file for SAM.
    if algorithm == "sam" and sam_qq_file:
        R('S <- x$S')
        jmath.R_fn(
            "bitmap", sam_qq_file, type="png256",
            height=1600, width=1600, units="px", res=300)
        jmath.R_fn("samr.plot", jmath.R_var("S"), sam_DELTA)
        jmath.R_fn("dev.off")

    # Convert this DataFrame into a Python object.  Columns of floats
    # can be StrVector objects if there are NA embedded within them.
    # NA are special objects of either type
    # rpy2.rinterface.NACharacterType or type
    # rpy2.rinterface.NARealType.
    tDATA_py = []
    header = [DATA_R.colnames[i] for i in range(DATA_R.ncol)]
    for zzz, col_R in enumerate(DATA_R):  # iterate over columns
        col_py = [col_R[i] for i in range(len(col_R))]

        if col_R.__class__.__name__ == "StrVector":
            pass
        elif col_R.__class__.__name__ == "FloatVector":
            col_py = [float(x) for x in col_py]
        elif col_R.__class__.__name__ == "IntVector":
            col_py = [int(x) for x in col_py]
        tDATA_py.append(col_py)
    DATA_py = jmath.transpose(tDATA_py)

    #handle = open('test01.txt', 'w')
    #for x in DATA_py:
    #    print >>handle, "\t".join(map(str, x))

    # Convert NA to None.
    for i in range(len(DATA_py)):
        for j in range(len(DATA_py[i])):
            if type(DATA_py[i][j]) in [
                rinterface.NACharacterType, rinterface.NARealType]:
                DATA_py[i][j] = None

    # Sort by increasing p-value, then decreasing fold change.
    name  = "p.value"
    direction = 1
    #if algorithm == "sam":
    #    name = "Score(d)"
    if name not in header:
        name = "Log_2 Fold Change"
        direction = -1
    assert name in header, 'I could not find the "%s" column.' % name

    I = header.index(name)
    #schwartz = [(direction*float(x[I]), x) for x in DATA_py]
    values = [x[I] for x in DATA_py]
    for i in range(len(values)):
        if values[i] is None:
            values[i] = direction*1E10
        else:
            values[i] = direction*float(values[i])
    schwartz = zip(values, DATA_py)
    schwartz.sort()
    DATA_py = [x[-1] for x in schwartz]

    # Filter based on user criteria.
    if fold_change is not None:
        log_2_fc = math.log(fold_change, 2)
        name = "Log_2 Fold Change"
        assert name in header, 'I could not find the "%s" column.' % name
        I = header.index(name)
        DATA_py = [x for x in DATA_py
                   if x[I] is not None and abs(x[I]) >= log_2_fc]
    if p_cutoff is not None:
        name  = "p.value"
        assert name in header, 'I could not find the "%s" column.' % name
        I = header.index(name)
        DATA_py = [x for x in DATA_py
                   if x[I] is not None and float(x[I]) < p_cutoff]
    if fdr_cutoff is not None:
        name  = "FDR"
        # This might be missing if all the genes have already been
        # filtered.
        #assert name in header, 'I could not find the "%s" column.' % name
        if name in header:
            I = header.index(name)
            DATA_py = [x for x in DATA_py
                       if x[I] is not None and float(x[I]) < fdr_cutoff]
    if bonf_cutoff is not None:
        name  = "Bonf"
        assert name in header, 'I could not find the "%s" column.' % name
        I = header.index(name)
        DATA_py = [x for x in DATA_py
                   if x[I] is not None and float(x[I]) < bonf_cutoff]

    ## If no significant genes, then don't produce any output.
    ##if not DATA_py:
    ##    return

    # Write to the outhandle.
    _write_matrix(outfile, header, DATA_py)
    # Don't close someone else's file handle.
    #outhandle.close()

    # Write out the gene sets in GMT format, if requested.
    if not gmt_file:
        return
    assert "Direction" in header, 'I could not find the "Direction" column.'
    assert "Gene ID" in header, 'I could not find the "Gene ID" column.'
    assert "Gene Name" in header, 'I could not find the "Gene Name" column.'
    I_direction = header.index("Direction")
    I_geneid = header.index("Gene ID")
    I_genename = header.index("Gene Name")

    # "Higher in <name1>"
    # "Higher in <name2>"
    possible_directions = ["Higher in %s" % name1, "Higher in %s" % name2]
    direction = [x[I_direction] for x in DATA_py]
    for x in direction:
        assert x.startswith("Higher in ")
        assert x in possible_directions
    samples = [x.replace("Higher in ", "") for x in direction]

    genesets = []  # list of (<SAMPLE>, [UP|DN])
    for s in samples:
        assert s in [name1, name2]
        # Make genesets relative to name2.  (Assume name1 is control).
        d = "UP"
        if s == name1:
            s, d = name2, "DN"
        genesets.append((s, d))
    genesets_all = sorted({}.fromkeys(genesets))

    outhandle = open(gmt_file, 'w')
    for geneset in genesets_all:
        sample, direct = geneset
        I = [i for (i, gs) in enumerate(genesets) if gs == geneset]
        gid = [DATA_py[i][I_geneid] for i in I]
        gn = [DATA_py[i][I_genename] for i in I]
        # gn might be float.  genesetlib expects array of strings.
        #import sys; sys.exit(0)
        gid = genesetlib.clean_genes(gid)
        gn = genesetlib.clean_genes(gn, delim=genename_delim)
        # <SAMPLE>_[ID|NAME]_[UP|DN]
        if gid:
            x = "%s_%s_%s" % (sample, "ID", direct)
            x = [x, "na"] +  gid
            print >>outhandle, "\t".join(x)
        if gn:
            x = "%s_%s_%s" % (sample, "NAME", direct)
            x = [x, "na"] + gn
            print >>outhandle, "\t".join(x)
    outhandle.close()


def _run_forked(
    outfile, gmt_file, algorithm, paired, MATRIX,
    geneid_header, genename_header, genename_delim,
    name1, name2, classes,
    fold_change, p_cutoff, fdr_cutoff, bonf_cutoff,
    sam_delta, sam_qq_file, edger_tagwise_dispersion, num_procs):
    import os
    import sys
    import tempfile

    tmpfile = pid = None
    try:
        x, tmpfile = tempfile.mkstemp(dir="."); os.close(x)
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)

        # Fork a subprocess, because some R libraries generate garbage
        # to the screen.
        r, w = os.pipe()
        pid = os.fork()
        if pid:   # Parent
            os.close(w)
            r = os.fdopen(r)
            for line in r:
                sys.stdout.write(line)   # output from R library
            sys.stdout.flush()
            os.waitpid(pid, 0)

            assert os.path.exists(tmpfile), "File not found: %s" % tmpfile
            if outfile is None:
                outfile = sys.stdout
            elif type(outfile) is type(""):
                outfile = open(outfile, 'w')
            for line in open(tmpfile):
                outfile.write(line)
            outfile.flush()
        else:     # Child
            os.close(r)
            w = os.fdopen(w, 'w')
            os.dup2(w.fileno(), sys.stdout.fileno())
            find_diffexp_genes(
                tmpfile, gmt_file, algorithm, paired, 
                MATRIX, geneid_header, genename_header, genename_delim,
                name1, name2, classes,
                fold_change, p_cutoff, fdr_cutoff, bonf_cutoff,
                sam_delta, sam_qq_file, edger_tagwise_dispersion, num_procs)
            sys.exit(0)
    finally:
        if pid:
            if tmpfile and os.path.exists(tmpfile):
                os.unlink(tmpfile)
    

def _run_not_forked(
    outfile, gmt_file, algorithm, paired, MATRIX,
    geneid_header, genename_header, genename_delim,
    name1, name2, classes,
    fold_change, p_cutoff, fdr_cutoff, bonf_cutoff,
    sam_delta, sam_qq_file, edger_tagwise_dispersion, num_procs):
    import sys

    outfile = outfile or sys.stdout
    find_diffexp_genes(
        outfile, gmt_file, algorithm, paired, 
        MATRIX, geneid_header, genename_header, genename_delim, 
        name1, name2, classes, fold_change, p_cutoff, fdr_cutoff, bonf_cutoff,
        sam_delta, sam_qq_file, edger_tagwise_dispersion, num_procs)
        

def main():
    import os
    import sys
    import argparse
    #from collections import Counter

    import arrayio
    from genomicode import arraysetlib
    #from genomicode import binreg
    #from genomicode import jmath
    
    parser = argparse.ArgumentParser(
        description="Find differentially expressed genes.")
    parser.add_argument("expression_file", help="Gene expression file.")

    # DESeq2 should have raw counts data.

    #parser.add_argument(
    #    "-l", "--log_the_data", 
    #    choices=["yes", "no", "auto"], default="auto",
    #    help="Log the data before analyzing.  "
    #    "Must be 'yes', 'no', or 'auto' (default).")
    #parser.add_argument(
    #    "--show_all_genes", default=False, action="store_true",
    #    help="List all the genes (default shows only p<0.05).")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of processors to use.")
    parser.add_argument(
        "--geneid_header", help="The column header of the gene IDs.")
    parser.add_argument(
        "--genename_header", help="The column header of the gene names.")
    parser.add_argument(
        "--num_header_cols", type=int,
        help="This number of columns are headers.  If not given, will guess.")
    parser.add_argument(
        "--gmt_file", help="Save the results in GMT format.")
    parser.add_argument(
        "--unfiltered_file", help="Save unfiltered results.")
    parser.add_argument(
        "--genename_delim",
        help="When writing gmt file, separate gene names based on this "
        "delimiter.")
    
    group = parser.add_argument_group(title="Class Labels")
    group.add_argument(
        "--cls_file", default=None, help="Class label file.")
    group.add_argument(
        "--indexes1", default=None,
        help="Which columns in class 1, e.g. 1-5,8 (1-based, inclusive).")
    group.add_argument(
        "--indexes2", default=None,
        help="Which columns in class 2, e.g. 1-5,8 (1-based, inclusive).  "
        "If not given, then will use all samples not in indexes1.")
    group.add_argument(
        "--indexes_include_headers", "--iih",
        default=False, action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    group.add_argument("--name1", default=None, help="Name for class 1.")
    group.add_argument("--name2", default=None, help="Name for class 2.")

    group = parser.add_argument_group(title="Cutoffs")
    group.add_argument(
        "--fold_change", type=float, default=None,
        help="Minimum change in gene expression (without log).  ")
    group.add_argument(
        "--p_cutoff", default=None, type=float,
        help="Only keep genes with p-value less than this value.")
    group.add_argument(
        "--fdr_cutoff", default=None, type=float,
        help="Only keep genes with FDR less than this value.")
    group.add_argument(
        "--bonf_cutoff", default=None, type=float,
        help="Only keep genes with BONFERRONI less than this value.")
    
    group = parser.add_argument_group(title="Algorithm-specific Parameters")
    group.add_argument(
        "--algorithm", 
        choices=["fold_change", "ttest", "sam", "ebayes", "deseq2", "edger"],
        default="ebayes", help="Which algorithm to use.")
    group.add_argument(
        "--paired", action="store_true",
        help="Do a paired analysis.  The same number of samples should be in "
        "each group.  Assumes that the samples should be paired in the same "
        "order that they occur in the file.  For example, if samples 1-5 are "
        "in group 1 and samples 6-10 are in group 2, then sample 1 will be "
        "paired with sample 6, 2 with 7, etc.")
    group.add_argument(
        "--sam_delta", type=float, default=1.0,
        help="DELTA parameter for SAM analysis.  (Default 1.0).")
    group.add_argument(
        "--sam_qq_file", help="File (PNG) for QQ plot for SAM.")
    group.add_argument(
        "--edger_tagwise_dispersion", action="store_true",
        help="Use tagwise rather than common dispersion.")


    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
    if args.num_header_cols is not None:
        assert args.num_header_cols > 0 and args.num_header_cols < 100
    assert args.num_procs >= 1 and args.num_procs < 100
    if args.fold_change is not None:
        assert args.fold_change >= 1E-10 and args.fold_change < 1000, \
               "Invalid fold change: %s" % args.fold_change
    if args.fdr_cutoff is not None:
        assert args.fdr_cutoff > 0.0 and args.fdr_cutoff <= 1.0
        assert args.p_cutoff is None, "Cannot have both FDR and p cutoff."
        assert args.bonf_cutoff is None, \
               "Cannot have both FDR and bonferroni cutoff."
    if args.p_cutoff is not None:
        assert args.bonf_cutoff is None, \
               "Cannot have both p and bonferroni cutoff."
        assert args.p_cutoff > 0.0 and args.p_cutoff <= 1.0
    if args.bonf_cutoff is not None:
        assert args.bonf_cutoff > 0.0 and args.bonf_cutoff <= 1.0
    if args.algorithm == "fold_change":
        assert not args.p_cutoff, \
               "Cannot use p-value cutoff with fold change algorithm."
        assert not args.fdr_cutoff, \
               "Cannot use fdr cutoff with fold change algorithm."
        assert not args.bonf_cutoff, \
               "Cannot use Bonferroni cutoff with fold change algorithm."
    
    assert args.sam_delta > 0 and args.sam_delta < 100
    if args.edger_tagwise_dispersion:
        assert args.algorithm == "edger", "tagwise dispersion only for edgeR"

    # Must have either the indexes or the cls_file, but not both.
    assert args.cls_file or args.indexes1, (
        "Must provide either CLS file or indexes.")
    if args.cls_file:
        assert os.path.exists(args.cls_file), \
               "File not found: %s" % args.cls_file
        assert not args.indexes1 and not args.indexes2
        assert not args.name1 and not args.name2
    # Don't need to check.  If name1 is missing, then a default will
    # be provided.
    #if args.indexes1:
    #    assert args.name1, "--name1 missing"
    #    assert args.name2, "--name2 missing"
    if args.indexes2:
        assert args.indexes1

    MATRIX = arrayio.read(args.expression_file, hcols=args.num_header_cols)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    #log_data = False
    #if args.log_the_data == "yes":
    #    log_data = True
    #elif args.log_the_data == "auto":
    #    log_data = not binreg.is_logged_array_data(MATRIX)
    #if log_data:
    #    MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)
    #    for i in range(len(MATRIX._X)):
    #        for j in range(len(MATRIX._X[i])):
    #            MATRIX._X[i][j] = max(MATRIX._X[i][j], 0)
    

    # Make a CLS file, if necessary.
    if args.cls_file:
        names, classes = arraysetlib.read_cls_file(args.cls_file)
        assert len(names) == 2, "I must have 2 classes."
        name1, name2 = names
        # Make sure classes has only 0, 1, or None.
        for i, x in enumerate(classes):
            if x == name1:
                classes[i] = 0
            elif x == name2:
                classes[i] = 1
            assert classes[i] in [None, 0, 1]
    else:
        x = arraysetlib.resolve_classes(
            MATRIX, args.indexes1, args.indexes2, args.indexes_include_headers,
            args.name1, args.name2)
        name1, name2, classes = x

    # Doesn't work.  Some algorithms need to do the filtering, because
    # it can change the multiple hypothesis correction.
    #args.unfiltered_file = None

    # Run the analysis
    #fn = _run_forked
    fn = _run_not_forked  # for debugging
    params = (
        sys.stdout, args.gmt_file, 
        args.algorithm, args.paired, MATRIX,
        args.geneid_header, args.genename_header, args.genename_delim, 
        name1, name2, classes,
        args.fold_change, args.p_cutoff, args.fdr_cutoff, args.bonf_cutoff,
        args.sam_delta, args.sam_qq_file, args.edger_tagwise_dispersion,
        args.num_procs)
    fn(*params)
    if args.unfiltered_file:
        params = (
            args.unfiltered_file, None, 
            args.algorithm, args.paired, MATRIX,
            args.geneid_header, args.genename_header, args.genename_delim, 
            name1, name2, classes,
            None, None, None, None, 
            args.sam_delta, args.sam_qq_file, args.edger_tagwise_dispersion,
            args.num_procs)
        fn(*params)
        
    
        

if __name__ == '__main__':
    main()

