#!/usr/bin/env python

import sys, os

# parse_geneset
# score_genes

def parse_geneset(geneset):
    # Return tuple of (positive_geneset, negative_geneset).
    # negative_geneset can be None.
    
    # geneset is in the format:
    # <positive_geneset>
    # <positive_geneset>,<negative_geneset>
    assert type(geneset) is type("")
    x = geneset.split(",")
    assert len(x) <= 2, "Unknown format for geneset: %s" % x
    pos, neg = x[0], None
    if len(x) == 2:
        neg = x[1]
    return pos, neg

def score_gene_set(MATRIX, matrix_name, name, pos_genes, neg_genes):
    from genomicode import matrixlib
    from genomicode import jmath

    all_genes = pos_genes + neg_genes
    x = matrixlib.find_best_row_header(MATRIX, all_genes)
    header, num_found, found, not_found = x
    assert num_found, "I could not find any genes in gene set %s." % name

    MATRIX_p = MATRIX.matrix(row=pos_genes, row_header=header)
    MATRIX_n = MATRIX.matrix(row=neg_genes, row_header=header)

    num_rows = MATRIX_p.nrow() + MATRIX_n.nrow()
    x = "Gene set %s contains %d genes and matched %d rows in %s."
    print x % (name, len(all_genes), num_rows, matrix_name)

    X_p = MATRIX_p._X
    X_n = []
    for x in MATRIX_n._X:
        x = [-x for x in x]
        X_n.append(x)
    X_pn = X_p + X_n
    score = jmath.mean(X_pn, byrow=False)

    return score

def main():
    import argparse
    import glob
    
    parser = argparse.ArgumentParser(
        description="Score a gene set on a gene expression data set.")

    parser.add_argument(
        "expression_files", nargs="+", help="Data set(s) to score.")

    # Assumes that there are no commas in names of gene sets.
    parser.add_argument(
        "--geneset_file", dest="geneset_files", action="append", default=[],
        help="File(s) with gene sets.  Should be in gmx or gmt format.")
    parser.add_argument(
        "-g", dest="gene_set", action="append", default=[],
        help="Name of the gene set to score.  If you want to score both "
        "the positively and negatively correlated genes, specify both "
        "gene sets using the format: <positive_geneset>,<negative_geneset>.  "
        "You can use this option multiple times to score more than one gene "
        "set.")
    parser.add_argument(
        "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_argument(
        "-o", dest="outfile", default=None, help="Name of file for results.")
    
    args = parser.parse_args()
    assert args.expression_files, \
           "Please specify an expression data set to score."
    expression_files = []
    for x in args.expression_files:
        x = glob.glob(x)
        expression_files.extend(x)
    for x in expression_files:
        assert os.path.exists(x), \
           "I could not find the expression file: %s" % x
    assert args.outfile, "Please specify the name of an outfile."
    assert args.geneset_files, "Please specify one or more geneset files."
    for x in args.geneset_files:
        assert os.path.exists(x), "I could not find the gene set file: %s" % x
    assert args.gene_set, "Please specify one or more gene sets to score."

    if args.libpath:
        sys.path = options.libpath + sys.path
    # Import after the library path is set.
    import time
    import arrayio
    from genomicode import genesetlib
    from genomicode import genepattern
    
    start_time = time.time()
    
    genepattern.fix_environ_path()

    msg = "Reading gene set file."
    if len(args.geneset_files) > 1:
        msg = "Reading gene set files."
    print msg; sys.stdout.flush()
    genesets = {}  # name -> list of genes
    for filename in args.geneset_files:
        for x in genesetlib.read_genesets(filename):
            name, description, genes = x
            assert name not in genesets, "Duplicate geneset: %s." % name
            genesets[name] = genes

    NAMES = []
    MATRICES = []
    for filename in expression_files:
        print "Reading gene expression file: %s." % filename
        sys.stdout.flush()
        x = os.path.split(filename)[1]
        NAMES.append(x)
        x = arrayio.read(filename)
        MATRICES.append(x)

    results = {}   # (matrix, geneset, sample) -> score
    for geneset in args.gene_set:
        pos_gs, neg_gs = parse_geneset(geneset)
        assert pos_gs in genesets, "I could not find gene set: %s" % pos_gs
        if neg_gs:
            assert neg_gs in genesets, "I could not find gene set: %s" % neg_gs
        gs_name = pos_gs
        if neg_gs:
            gs_name = "%s/%s" % (pos_gs, neg_gs)
            
        pos_genes = genesets[pos_gs]
        neg_genes = genesets.get(neg_gs, [])

        for matrix_name, MATRIX in zip(NAMES, MATRICES):
            scores = score_gene_set(
                MATRIX, matrix_name, gs_name, pos_genes, neg_genes)
            sample_names = MATRIX.col_names(arrayio.COL_ID)
            assert len(sample_names) == len(scores)
            for (sample, score) in zip(sample_names, scores):
                key = matrix_name, gs_name, sample
                assert key not in results, "Duplicate: %s" % key
                results[key] = score

    x = [(x[0], x[2]) for x in results]
    x = sorted({}.fromkeys(x))
    all_matrix_samples = x
    x = [x[1] for x in results]
    x = sorted({}.fromkeys(x))
    all_genesets = x

    outhandle = open(args.outfile, 'w')
    header = ["FILE", "SAMPLE"] + all_genesets
    print >>outhandle, "\t".join(header)
    for x in all_matrix_samples:
        matrix, sample = x
        scores = [results[(matrix, x, sample)] for x in all_genesets]
        x = [matrix, sample] + scores
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
    
    print "Done."
    
if __name__ == '__main__':
    main()
