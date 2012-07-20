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

def score_gene_set(MATRIX, name, pos_genes, neg_genes):
    from genomicode import matrixlib
    from genomicode import jmath

    all_genes = pos_genes + neg_genes
    x = matrixlib.find_best_row_header(MATRIX, all_genes)
    header, num_found, found, not_found = x
    assert num_found, "I could not find any genes in gene set %s." % name

    MATRIX_p = MATRIX.matrix(row=pos_genes, row_header=header)
    MATRIX_n = MATRIX.matrix(row=neg_genes, row_header=header)

    num_rows = MATRIX_p.nrow() + MATRIX_n.nrow()
    x = "Gene set %s contains %d genes and matched %d rows in " + \
        "the expression data."
    print x % (name, len(all_genes), num_rows)

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
    
    parser = argparse.ArgumentParser(
        description="Score a gene set on a gene expression data set.")

    parser.add_argument("geneset_file", help="File with gene sets.  "
                        "Should be in gmx or gmt format.")
    parser.add_argument("expression_file", help="Data set to score.")

    # Assumes that there are no commas in names of gene sets.
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
    assert os.path.exists(args.geneset_file), \
           "I could not find the gene set file: %s" % args.geneset_file
    assert os.path.exists(args.expression_file), \
           "I could not find the expression file: %s" % args.expression_file
    assert args.outfile, "Please specify the name of an outfile."

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

    print "Reading gene set file."; sys.stdout.flush()
    genesets = {}  # name -> list of genes
    for x in genesetlib.read_genesets(args.geneset_file):
        name, description, genes = x
        genesets[name] = genes
    
    print "Reading gene expression file."; sys.stdout.flush()
    MATRIX = arrayio.read(args.expression_file)

    results = []   # list of (geneset, list of scores)
    for geneset in args.gene_set:
        pos_gs, neg_gs = parse_geneset(geneset)
        assert pos_gs in genesets, "I could not find gene set: %s" % pos_gs
        if neg_gs:
            assert neg_gs in genesets, "I could not find gene set: %s" % neg_gs
        name = pos_gs
        if neg_gs:
            name = "%s/%s" % (pos_gs, neg_gs)
            
        pos_genes = genesets[pos_gs]
        neg_genes = genesets.get(neg_gs, [])
        scores = score_gene_set(MATRIX, name, pos_genes, neg_genes)
        results.append((name, scores))

    outhandle = open(args.outfile, 'w')
    sample_names = MATRIX.col_names(arrayio.COL_ID)
    header = ["SAMPLE"] + [x[0] for x in results]
    print >>outhandle, "\t".join(header)
    for i in range(len(sample_names)):
        x = [sample_names[i]]
        for (name, scores) in results:
            assert len(scores) == len(sample_names), "%d %d" % (
                len(scores), len(sample_names))
            x.append(scores[i])
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
    
    print "Done."
    
if __name__ == '__main__':
    main()
