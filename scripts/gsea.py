#!/usr/bin/env python

# Still some problems with the platforms.
# 1.  Convert to GCT may discard the correct annotations.
# 2.  What if the file already contains gene symbols?  I'm not sure
# how to specify this for GenePattern.  For the website, we may be
# able to leave it blank.  But this doesn't seem to work through the R
# interface.



# Mapping from arrayplatformlib to GenePattern.
platform2gpplatform = {
    "HG_U133A" : "HG_U133A.chip",
    "HG_U133A_2" : "HG_U133A_2.chip",
    "HG_U133_Plus_2" : "HG_U133_Plus_2.chip",
    #"Entrez_symbol_human" : None,
    # Missing entrez_ID_human
    }

def guess_chip_platform(M):
    # Return the GenePattern chip.platform for this matrix.
    from genomicode import arrayplatformlib

    platform = arrayplatformlib.identify_platform_of_matrix(M)
    if platform is None:
        return None
    assert platform in platform2gpplatform, \
        "I don't know how to convert %s to a GenePattern platform." % platform
    chipname = platform2gpplatform.get(platform)
    return chipname

DATABASE2GENESET = {
    "positional" : "c1.all.v3.0.symbols.gmt",
    "curated" : "c2.all.v3.0.symbols.gmt",
    "curated:canonical" : "c2.cp.v3.0.symbols.gmt",
    "curated:biocarta" : "c2.cp.biocarta.v3.0.symbols.gmt",
    "curated:kegg" : "c2.cp.kegg.v3.0.symbols.gmt",
    "curated:reactome" : "c2.cp.reactome.v3.0.symbols.gmt",
    "motif" : "c3.all.v3.0.symbols.gmt",
    "motif:tfactor" : "c3.tft.v3.0.symbols.gmt",
    "computational" : "c4.all.v3.0.symbols.gmt",
    "gene_ontology" : "c5.all.v3.0.symbols.gmt",
    "gene_ontology:process" : "c5.bp.v3.0.symbols.gmt",
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
        "--indexes_include_headers", default=False, action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    group.add_argument("--name1", default=None, help="Name for group 1.")
    group.add_argument("--name2", default=None, help="Name for group 2.")

    group = parser.add_argument_group(title="Other Parameters")
    group.add_argument(
        "--platform", default=None,
        help="The platform (GenePattern chip) of the expression data, "
        "e.g. HG_U133A_2.chip")
    x = sorted(DATABASE2GENESET)
    x = [x.replace(DEFAULT_DATABASE, "%s (DEFAULT)" % x) for x in x]
    x = ", ".join(x)
    x = "Which database to search.  Possible values are: %s." % x
    group.add_argument("--database", default=DEFAULT_DATABASE, help=x)
    group.add_argument("--database_file", default=None,
        help="Search a GMT or GMX file instead of the default databases.")
    # GenePattern uses the platform (chip file) to convert the IDs in
    # the Gene Symbols.  If the IDs are already gene symbols, then an
    # arbitrary platform should be given (since it requires this
    # parameter), and collapse should be turned off.
    group.add_argument(
        "--no_collapse_dataset", default=False, action="store_true",
        help="Do not 1) convert gene IDs to gene symbols, and do not 2) "
        "collapse duplicate gene symbols.  Set this if the gene IDs are "
        "already unique gene symbols.  Also, can use this if you "
        "provide the database_file and the gene IDs match the ones "
        "in our gene expression file.")
    
    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
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

    # BUG: Should check the names of the samples to make sure they
    # don't contain spaces or other punctuation.  GSEA seems to be
    # sensitive to these things.

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
    if platform is None and not database_file:
        platform = guess_chip_platform(MATRIX)
        assert platform, "I could not find a platform for this file."

    # Set up file names.
    opj = os.path.join
    x = os.path.splitext(os.path.split(args.expression_file)[1])[0]
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
        }
    # platform is required, even if collapse.dataset is false.  If no
    # platform is given, then specify a default one.
    params["chip.platform"] = platform
    if params["chip.platform"] is None:
        params["chip.platform"] = "HG_U133A.chip"
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
    #for x in cmd:
    #    print x
    #import sys; sys.exit(0)

    # Run the analysis in the outpath.  GSEA leaves a file
    # "System.out" in the current directory.
    cwd = os.getcwd()
    try:
        os.chdir(args.outpath)
        p = subprocess.Popen(
            cmd, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            close_fds=True)
        p.wait()
    finally:
        os.chdir(cwd)

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
