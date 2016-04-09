#!/usr/bin/env python

# Possibilities:
# 1.  Search database on GSEA.  Use some sort of gene ID.
#     A.  If provide platform, then use it.
#     B.  Try to guess platform.
#     C.  If Gene Symbol, turn off collapse_dataset.
#     D.  If Gene ID, convert to Gene Symbol (not implemented).
# 2.  Search own database (with database_file).
#     A.  collapse_dataset should be off.
#     B.  IDs in expression should match IDs in database.

# GenePattern uses the platform (chip file) to convert the IDs in the
# Gene Symbols.  If the IDs are already gene symbols, then no platform
# should be given (since it requires this parameter) and collapse
# should be turned off.


# Still some problems with the platforms.
# 1.  Convert to GCT may discard the correct annotations.


# Mapping from arrayplatformlib to GenePattern.
platform2gpplatform = {
    "HG_U133A" : "HG_U133A.chip",
    "HG_U133A_2" : "HG_U133A_2.chip",
    "HG_U133_Plus_2" : "HG_U133_Plus_2.chip",
    "Entrez_symbol_human" : None,
    "Entrez_Symbol_human" : None,
    "entrez_ID_symbol_human" : None,
    # Missing entrez_ID_human
    }


def guess_chip_platform(M, min_match_score):
    # Return the GenePattern chip.platform for this matrix.
    from genomicode import arrayplatformlib as apl

    #platform = arrayplatformlib.identify_platform_of_matrix(M)
    #assert platform, "I could not guess the platform for this file."
    x = apl.score_matrix(M, min_score=0.01)
    assert x, "I could not guess the platform for this file.%s" % x
    best_score = x[0]
    
    x = ""
    if best_score.max_score > 0:
        x = "  The closest was %s (%.2f)." % (
            best_score.platform_name, best_score.max_score)
    assert best_score.max_score > min_match_score, \
           "I could not guess the platform for this file.%s" % x
    assert best_score.platform_name in platform2gpplatform, \
           "I don't know how to convert %s to a GenePattern platform." % \
           best_score.platform_name
    chipname = platform2gpplatform.get(best_score.platform_name)
    return chipname

DATABASE2GENESET = {
    ## "positional" : "c1.all.v3.0.symbols.gmt",
    ## "curated" : "c2.all.v3.0.symbols.gmt",
    ## "curated:canonical" : "c2.cp.v3.0.symbols.gmt",
    ## "curated:biocarta" : "c2.cp.biocarta.v3.0.symbols.gmt",
    ## "curated:kegg" : "c2.cp.kegg.v3.0.symbols.gmt",
    ## "curated:reactome" : "c2.cp.reactome.v3.0.symbols.gmt",
    ## "motif" : "c3.all.v3.0.symbols.gmt",
    ## "motif:tfactor" : "c3.tft.v3.0.symbols.gmt",
    ## "computational" : "c4.all.v3.0.symbols.gmt",
    ## "gene_ontology" : "c5.all.v3.0.symbols.gmt",
    ## "gene_ontology:process" : "c5.bp.v3.0.symbols.gmt",
    
    ## "positional" : "c1.all.v4.0.symbols.gmt",
    ## "curated" : "c2.all.v5.0.symbols.gmt",
    ## "curated:canonical" : "c2.cp.v4.0.symbols.gmt",
    ## "curated:biocarta" : "c2.cp.biocarta.v4.0.symbols.gmt",
    ## "curated:kegg" : "c2.cp.kegg.v4.0.symbols.gmt",
    ## "curated:reactome" : "c2.cp.reactome.v4.0.symbols.gmt",
    ## "motif" : "c3.all.v4.0.symbols.gmt",
    ## "motif:tfactor" : "c3.tft.v4.0.symbols.gmt",
    ## "computational" : "c4.all.v4.0.symbols.gmt",
    ## "gene_ontology" : "c5.all.v4.0.symbols.gmt",
    ## "gene_ontology:process" : "c5.bp.v4.0.symbols.gmt",
    
    "positional" : "c1.all.v5.0.symbols.gmt",
    "curated" : "c2.all.v5.0.symbols.gmt",
    "curated:canonical" : "c2.cp.v5.0.symbols.gmt",
    "curated:biocarta" : "c2.cp.biocarta.v5.0.symbols.gmt",
    "curated:kegg" : "c2.cp.kegg.v5.0.symbols.gmt",
    "curated:reactome" : "c2.cp.reactome.v5.0.symbols.gmt",
    "motif" : "c3.all.v5.0.symbols.gmt",
    "motif:tfactor" : "c3.tft.v5.0.symbols.gmt",
    "computational" : "c4.all.v5.0.symbols.gmt",
    "gene_ontology" : "c5.all.v5.0.symbols.gmt",
    "gene_ontology:process" : "c5.bp.v5.0.symbols.gmt",
    }
DEFAULT_DATABASE = "gene_ontology:process"

def format_gene_set_database(database):
    # Valid values for GenePattern gene.set.database:
    # c1.all.v3.0.symbols.gmt [Positional]
    # c2.all.v3.0.symbols.gmt [Curated]
    # c2.cgp.v3.0.symbols.gmt [Curated]        chemical & genetic pertubations
    # c2.cp.v3.0.symbols.gmt [Curated]            canonical pathways
    # c2.cp.biocarta.v3.0.symbols.gmt [Curated]   Biocarta
    # c2.cp.kegg.v3.0.symbols.gmt [Curated]       KEGG
    # c2.cp.reactome.v3.0.symbols.gmt [Curated]   Reactome
    # c3.all.v3.0.symbols.gmt [Motif]
    # c3.mir.v3.0.symbols.gmt [Motif]             miRNA targets
    # c3.tft.v3.0.symbols.gmt [Motif]             transcription factor targets
    # c4.all.v3.0.symbols.gmt [Computational]
    # c4.cgn.v3.0.symbols.gmt [Computational]     cancer gene neighborhoods
    # c4.cm.v3.0.symbols.gmt [Computational]      cancer modules
    # c5.all.v3.0.symbols.gmt [Gene Ontology]
    # c5.bp.v3.0.symbols.gmt [Gene Ontology]      biological process
    # c5.cc.v3.0.symbols.gmt [Gene Ontology]      cellular component
    # c5.mf.v3.0.symbols.gmt [Gene Ontology]      molecular function
    assert database in DATABASE2GENESET, "Unknown database: %s" % database
    return DATABASE2GENESET[database]


def fix_class_order(MATRIX, name1, name2, classes):
    # Make sure classes are in the right order.  If not, reorder the
    # matrix so that the samples for the first class come first.

    assert classes

    # Convert classes to 0/1.
    clean = []
    for c in classes:
        if c in ["0", "1"]:
            c = int(c)
        elif c in [name1, name2]:
            c = int(c == name2)
        assert c in [0, 1], "Unknown class: %s" % c
        clean.append(c)
    classes = clean

    # If the classes are in the right order, don't need to do anything.
    if classes[0] == 0:
        x = MATRIX, name1, name2, classes
        return x
    
    indexes0, indexes1 = [], []
    for i, c in enumerate(classes):
        assert c in [0, 1]
        if c == 0:
            indexes0.append(i)
        else:
            indexes1.append(i)
    O = indexes0 + indexes1

    classes = [classes[i] for i in O]
    MATRIX = MATRIX.matrix(None, O)

    x = MATRIX, name1, name2, classes
    return x


def check_matrix(X):
    import re
    import arrayio
    #from genomicode import hashlib

    assert arrayio.gct_format.is_matrix(X)

    # Make sure gene IDs (NAME) is unique and non-empty.
    assert X.row_names()[0].upper() == "NAME", \
           "Header of first column should be: NAME"
    seen = {}
    for i, name in enumerate(X.row_names("NAME")):
        assert name.strip(), "Empty gene ID in row %d." % (i+1)
        assert name not in seen, "Duplicate gene ID: %s" % name
        seen[name] = 1

    # Make sure sample names don't contain spaces or other
    # punctuation.  GSEA seems to be sensitive to these things.
    sample_names = X.col_names(arrayio.tdf.SAMPLE_NAME)
    bad_names = []
    for i, name in enumerate(sample_names):
        if not name:
            bad_names.append("<blank>")
        elif re.search("[^a-zA-Z0-9_-]", name):
            bad_names.append(name)
            
    # If there are bad names, try to fix them.
    #if bad_names:
    #    sample_names_h = [hashlib.hash_var(x) for x in sample_names]
    #    # If there are no duplicates, use these sample names.
    #    raise NotImplementedError
    assert not bad_names, "Bad sample name: %s" % ", ".join(bad_names)

    # Make sure sample names are unique.
    seen = {}
    for i, name in enumerate(sample_names):
        assert name not in seen, "Duplicate sample name: %s" % name
        seen[name] = 1


def check_classes(classes, permutation_type):
    assert classes

    class2count = {}
    for x in classes:
        if x not in class2count:
            class2count[x] = 0
        class2count[x] += 1
    
    # Make sure there are exactly 2 classes.
    assert len(class2count) != 1, "Only 1 class given."
    assert len(class2count) == 2, "Cannot have more than 2 classes."

    # If permutation_type=phenotype, need at least 3 samples per class.
    if permutation_type == "phenotype":
        counts = class2count.values()
        assert min(counts) >= 3, \
               "Need at least 3 samples for phenotype permutations."
        

def main():
    import os
    import argparse
    import subprocess
    import StringIO
    import zipfile
    import shutil

    import arrayio
    from genomicode import config
    from genomicode import arraysetlib
    
    parser = argparse.ArgumentParser(description="Do a GSEA analysis.")
    parser.add_argument("expression_file", help="Gene expression file.")
    parser.add_argument("outpath", help="Where to save the files.")
    parser.add_argument(
        "--clobber", default=False, action="store_true",
        help="Overwrite outpath, if it already exists.")
    parser.add_argument(
        "--dry_run", default=False, action="store_true",
        help="Set up the file, but do not run GSEA.")
    
    group = parser.add_argument_group(title="Class Labels")
    group.add_argument(
        "--cls_file", default=None, help="Class label file.")
    group.add_argument(
        "--indexes1", default=None,
        help="Which columns in group 1, E.g. 1-5,8 (1-based, "
        "inclusive).  All other samples will be in group 2.")
    group.add_argument(
        "--indexes_include_headers", "--iih", action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    group.add_argument("--name1", default=None, help="Name for group 1.")
    group.add_argument("--name2", default=None, help="Name for group 2.")

    group = parser.add_argument_group(title="Other Parameters")
    group.add_argument(
        "--platform", default=None,
        help="The platform (GenePattern chip) of the expression data, "
        "e.g. HG_U133A_2.chip.  You should leave this blank if the IDs "
        "in the gene expression data set are gene symbols.  Allowed "
        "values can be found on the GenePattern/GSEA website.")
    group.add_argument(
        '--min_match_score', default=0.95, type=float,
        help="When trying to identify the rows of a matrix or geneset, "
        "require at least this portion of the IDs to be recognized.")
    x = sorted(DATABASE2GENESET)
    x = [x.replace(DEFAULT_DATABASE, "%s (DEFAULT)" % x) for x in x]
    x = ", ".join(x)
    x = "Which database to search.  Possible values are: %s." % x
    group.add_argument("--database", default=DEFAULT_DATABASE, help=x)
    group.add_argument("--database_file", default=None,
        help="Search a GMT or GMX file instead of the default databases.")
    group.add_argument(
        "--no_collapse_dataset", default=False, action="store_true",
        help="Do not 1) convert gene IDs to gene symbols, and do not 2) "
        "collapse duplicate gene symbols.  Set this if the gene IDs are "
        "already unique gene symbols.  Also, can use this if you "
        "provide the database_file and the gene IDs match the ones "
        "in our gene expression file.")

    # phenotype is more accurate.  But if only 2 samples, need to be
    # gene_set.  (Not sure about 3 samples?  Where is the limit?)
    group.add_argument(
        "--permutation_type", default="phenotype",
        choices=["phenotype", "gene_set"],
        help="Default is phenotype.  With <= 6 samples, recommend using "
        "gene_set instead.")
    
    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
    
    assert type(args.min_match_score) is type(0.0)
    assert args.min_match_score > 0.2, "min_match_score too low"
    assert args.min_match_score <= 1.0, "min_match_score too high"

    # Must have either the indexes or the cls_file, but not both.
    assert args.cls_file or args.indexes1, (
        "Must provide either CLS file or the indexes for one group.")
    assert not (args.cls_file and args.indexes1), (
        "Cannot provide both a CLS file and the indexes.")
    if args.cls_file:
        assert os.path.exists(args.cls_file), \
            "File not found: %s" % args.cls_file
        assert not args.name1
        assert not args.name2
    assert args.outpath, "Please specify an outpath."
    assert not os.path.exists(args.outpath) or args.clobber, \
        "Outpath %s already exists." % args.outpath
    if os.path.exists(args.outpath):
        shutil.rmtree(args.outpath)
    os.mkdir(args.outpath)

    MATRIX = arrayio.read(args.expression_file)

    # Make a CLS file, if necessary.
    if args.cls_file:
        names, classes = arraysetlib.read_cls_file(args.cls_file)
        assert len(names) == 2, "I must have 2 classes."
        name1, name2 = names
    else:
        x = arraysetlib.resolve_classes(
            MATRIX, args.indexes1, None, args.indexes_include_headers,
            args.name1, args.name2)
        name1, name2, classes = x

    x = fix_class_order(MATRIX, name1, name2, classes)
    MATRIX, name1, name2, classes = x

    handle = StringIO.StringIO()
    arraysetlib.write_cls_file(handle, name1, name2, classes)
    cls_data = handle.getvalue()

    # Convert the format after making CLS file, or else args.indexes1
    # with args.indexes_include_headers might be off.
    # BUG: What if the conversion to GCT discards the proper platform
    # of this matrix?
    MATRIX = arrayio.convert(MATRIX, to_format=arrayio.gct_format)

    database_file = None
    if args.database_file:
        database_file = os.path.realpath(args.database_file)
        assert os.path.exists(database_file), (
            "I could not find file: %s" % database_file)
    # Required, even if gene.sets.database.file is given.
    gene_set_database = format_gene_set_database(args.database)
    

    # If no database file is given, then we need to know the platform
    # for the expression data.  (If one is given, we let the user make
    # sure the platforms match.)
    platform = args.platform
    if platform is None and not database_file and not args.no_collapse_dataset:
        platform = guess_chip_platform(MATRIX, args.min_match_score)
        # If gene symbols already provided, then turn off collapse_dataset.
        if platform is None:  # Gene Symbol
            args.no_collapse_dataset = True

    # Do some sanity checking to make sure imputs are reasonable.
    check_matrix(MATRIX)
    check_classes(classes, args.permutation_type)
    
    # Set up file names.
    opj = os.path.join
    x = os.path.split(args.expression_file)[1]
    if x.lower().endswith(".gz"):
        x = x[:-3]
    x = os.path.splitext(x)[0]
    x = x.replace(" ", "_")  # GenePattern cannot work with spaces.
    assert x, "empty file name"
    gct_file = "%s.gct" % x
    cls_file = "%s.cls" % x
    out_file = "%s.zip" % x  # GenePattern saves output to <file>.zip.
    gct_full = opj(args.outpath, gct_file)
    cls_full = opj(args.outpath, cls_file)
    out_full = opj(args.outpath, out_file)
    if database_file:
        db_file = os.path.split(database_file)[1]
        db_full = opj(args.outpath, db_file)

    # Write the gene expression, class label, and database files.  It
    # is better to have local copies of the files.  It is unclear how
    # to upload files to GenePattern if the file names have spaces in
    # them.  Get around this by making all the files local.
    arrayio.gct_format.write(MATRIX, open(gct_full, 'w'))
    open(cls_full, 'w').write(cls_data)
    if database_file:
        open(db_full, 'w').write(open(database_file).read())

    collapse_dataset = "true"
    if args.no_collapse_dataset:
        collapse_dataset = "false"

    # Set up the analysis.
    params = {
        "expression.dataset" : gct_file,
        "phenotype.labels" : cls_file,
        "collapse.dataset" : collapse_dataset,
        "permutation.type" : args.permutation_type,
        }
    # platform is required, even if collapse.dataset is false.  If no
    # platform is given, then specify a default one.
    #CHIP_PLATFORM = "chip.platform"
    CHIP_PLATFORM = "chip.platform.file"
    params[CHIP_PLATFORM] = platform
    if params[CHIP_PLATFORM] is None:
        params[CHIP_PLATFORM] = "HG_U133A.chip"
    if database_file:
        params["gene.sets.database.file"] = db_file
    # Required, even if gene.sets.database.file is given.
    params["gene.sets.database"] = gene_set_database

    if args.dry_run:
        return
 
    cmd = [
        config.genepattern,
        "-o", ".",
        "GSEA"
        ]
    for (key, value) in reversed(list(params.iteritems())):
        x = ["--parameters", "%s:%s" % (key, value)]
        cmd.extend(x)
    #print " ".join(cmd)
    #import sys; sys.exit(0)

    # Run the analysis in the outpath.  GSEA leaves a file
    # "System.out" in the current directory.
    cwd = os.getcwd()
    try:
        os.chdir(args.outpath)
        p = subprocess.Popen(
            cmd, bufsize=0, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        # Check for errors in the output.
        x = r.read()
        # Get rid of GenePattern garbage.
        data = x.replace("Loading required package: rJava", "")
        p.wait()
    finally:
        os.chdir(cwd)
    x = data.strip()
    # rpy2 generates UserWarnings for some reason.
    # /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/
    #   python2.7/site-packages/rpy2/robjects/functions.py:106: UserWarning:
    #       res = super(Function, self).__call__(*new_args, **new_kwargs)
    if x.find("UserWarning") >= 0 and x.endswith("(*new_args, **new_kwargs)"):
        # Ignore this UserWarning.
        x = ""
    assert not x, "%s\n%s" % (cmd, data)


    error_file = os.path.join(args.outpath, "stderr.txt")
    assert not os.path.exists(error_file), (
        "Error generated by GenePattern:\n%s" % open(error_file).read())

    # Unzip the zipped results in the outpath.
    assert os.path.exists(out_full), "ERROR: Output file is missing [%s]." % \
        out_full
    zfile = zipfile.ZipFile(out_full)
    zfile.extractall(args.outpath)
    os.unlink(out_full)
                
    x = os.path.join(args.outpath, "index.html")
    assert os.path.exists(x), "I could not find the GSEA output: %s" % x
                                    
if __name__ == '__main__':
    main()
