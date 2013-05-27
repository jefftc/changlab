#!/usr/bin/env python


def choose_gene_names(MATRIX):
    # Return tuple of (header for gene_id, header for gene_names).
    # Either of the headers can be None.
    from genomicode import arrayplatformlib

    if not MATRIX.row_names():
        return None, None
    if len(MATRIX.row_names()) == 1:
        return MATRIX.row_names()[0], None
    if len(MATRIX.row_names()) == 2:
        return MATRIX.row_names()

    # By default, set as the first two rows.
    geneid_header, genename_header = MATRIX.row_names()[:2]

    # If I can find a better header for the annotation, then use it.
    x = arrayplatformlib.score_all_platforms_of_matrix(MATRIX)
    for header, platform, score in x:
        if score < 0.5:
            continue
        if platform == "entrez_ID_symbol_human":
            genename_header = header
            
    return geneid_header, genename_header
    
    
def find_diffexp_genes(
    outfile, gmt_file, algorithm, MATRIX, name1, name2, classes, fold_change,
    DELTA, num_procs):
    # classes must be 0, 1, None.
    import os
    
    from genomicode import config
    from genomicode import jmath
    from genomicode import genesetlib

    algorithm2function = {
        "ttest" : "find.de.genes.ttest",
        "sam" : "find.de.genes.sam",
        "ebayes" : "find.de.genes.ebayes",
        }
    assert algorithm in algorithm2function, "Unknown algorithm: %s" % algorithm

    # Select the relevant columns from MATRIX.
    I = [i for (i, x) in enumerate(classes) if x in [0, 1]]
    assert len(I)
    MATRIX = MATRIX.matrix(None, I)
    classes = [classes[i] for i in I]
    
    # Make sure at least 2 of each class.
    counts = {}
    for x in classes:
        counts[x] = counts.get(x, 0) + 1
    assert sorted(counts) == [0, 1], "Only one class represented."
    assert counts[0] >= 2, "There must be at least 2 of each class."
    assert counts[1] >= 2, "There must be at least 2 of each class."

    names = [name1, name2]
    X = MATRIX._X
    Y = [names[x] for x in classes]
    sample_name = None
    if MATRIX.col_names():
        sample_name = MATRIX.col_names(MATRIX.col_names()[0])

    x = choose_gene_names(MATRIX)
    geneid_header, genename_header = x

    R = jmath.start_R()
    de_lib = os.path.join(config.changlab_Rlib, "diffexp.R")
    stat_lib = os.path.join(config.changlab_Rlib, "statlib.R")
    assert os.path.exists(de_lib)
    assert os.path.exists(stat_lib)
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
        args.append("%g" % DELTA)
    if geneid:
        args.append("geneid=geneid")
    if genenames:
        args.append("genenames=genenames")
    if fold_change is not None:
        args.append("FOLD.CHANGE=%g" % fold_change)
    if algorithm == "ttest":
        args.append("NPROCS=%d" % num_procs)  # t-test only

    fn = algorithm2function[algorithm]
    x = ", ".join(args)
    R("x <- %s(%s)" % (fn, x))
    R("DATA <- x$DATA")
    DATA_R = R["DATA"]

    # Convert this DataFrame into a Python object.
    tDATA_py = []
    for col_R in DATA_R:  # iterate over columns
        col_py = [col_R[x] for x in range(len(col_R))]
        # What if there is NA?
        if col_R.__class__.__name__ == "StrVector":
            pass
        elif col_R.__class__.__name__ == "FloatVector":
            col_py = [float(x) for x in col_py]
        elif col_R.__class__.__name__ == "IntVector":
            col_py = [int(x) for x in col_py]
        tDATA_py.append(col_py)
    DATA_py = jmath.transpose(tDATA_py)
    header = [DATA_R.colnames[x] for x in range(DATA_R.ncol)]

    # Write to the outhandle.
    outhandle = open(outfile, 'w')
    print >>outhandle, "\t".join(header)
    outhandle.flush()
    for x in DATA_py:
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
        outhandle.flush()
    outhandle.close()

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
        gid = genesetlib.clean_genes(gid)
        gn = genesetlib.clean_genes(gn)
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
    

def main():
    import os
    import sys
    import argparse
    import tempfile

    import arrayio
    from genomicode import arraysetlib
    from genomicode import binreg
    from genomicode import jmath
    
    parser = argparse.ArgumentParser(
        description="Find differentially expressed genes.")
    parser.add_argument("expression_file", help="Gene expression file.")

    parser.add_argument(
        "-l", "--log_the_data", dest="log_the_data", 
        choices=["yes", "no", "auto"], default="auto",
        help="Log the data before analyzing.  "
        "Must be 'yes', 'no', or 'auto' (default).")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of processors to use.")
    parser.add_argument(
        "--gmt_file", help="Save the results in GMT format.")
    
    group = parser.add_argument_group(title="Algorithm Parameters")
    group.add_argument(
        "--algorithm", dest="algorithm", 
        choices=["ttest", "sam", "ebayes"], default="ebayes",
        help="Which algorithm to use.")
    group.add_argument(
        "--fold_change", type=float, default=None,
        help="Minimum change in gene expression.")
    group.add_argument(
        "--DELTA", type=float, default=1.0,
        help="DELTA parameter for SAM analysis.  (Default 1.0).")

    group = parser.add_argument_group(title="Class Labels")
    group.add_argument(
        "--cls_file", default=None, help="Class label file.")
    group.add_argument(
        "--indexes1", default=None,
        help="Which columns in class 1, E.g. 1-5,8 (1-based, inclusive).")
    group.add_argument(
        "--indexes2", default=None,
        help="Which columns in class 2, E.g. 1-5,8 (1-based, inclusive).  "
        "If not given, then will use all samples not in indexes1.")
    group.add_argument(
        "--indexes_include_headers", default=False, action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    group.add_argument("--name1", default=None, help="Name for class 1.")
    group.add_argument("--name2", default=None, help="Name for class 2.")

    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
    assert args.num_procs >= 1 and args.num_procs < 100
    if args.fold_change is not None:
        assert args.fold_change >= 0 and args.fold_change < 1000
    assert args.DELTA > 0 and args.DELTA < 100
    
    # Must have either the indexes or the cls_file, but not both.
    assert args.cls_file or args.indexes1, (
        "Must provide either CLS file or indexes.")
    if args.cls_file:
        assert os.path.exists(args.cls_file), \
               "File not found: %s" % args.cls_file
        assert not args.indexes1 and not args.indexes2
        assert not args.name1 and not args.name2
    if args.indexes1:
        assert args.name1 and args.name2
    if args.indexes2:
        assert args.indexes1

    MATRIX = arrayio.read(args.expression_file)

    log_data = False
    if args.log_the_data == "yes":
        log_data = True
    elif args.log_the_data == "auto":
        log_data = not binreg.is_logged_array_data(MATRIX)
    if log_data:
        MATRIX._X = jmath.log(MATRIX._X, base=2, safe=1)
    

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

    # Run the analysis.
    outfile = None
    try:
        x, outfile = tempfile.mkstemp(dir="."); os.close(x)

        # Fork a subprocess, because some R libraries generate garbage
        # to the screen.
        r, w = os.pipe()
        pid = os.fork()
        if pid:   # Parent
            os.close(w)
            r = os.fdopen(r)
            for line in r:
                sys.stdout.write(line)   # output from R library
                pass
            os.waitpid(pid, 0)

            assert os.path.exists(outfile), "failed"
            for line in open(outfile):
                sys.stdout.write(line)
        else:     # Child
            os.close(r)
            w = os.fdopen(w, 'w')
            os.dup2(w.fileno(), sys.stdout.fileno())
            find_diffexp_genes(
                outfile, args.gmt_file, args.algorithm, MATRIX,
                name1, name2, classes, args.fold_change, args.DELTA,
                args.num_procs)
            sys.exit(0)
    finally:
        if pid:
            if outfile and os.path.exists(outfile):
                os.unlink(outfile)
        

if __name__ == '__main__':
    main()
