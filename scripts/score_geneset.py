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

def _score_gene_set_h(MATRIX, matrix_name, name, pos_genes, neg_genes, lock):
    from genomicode import matrixlib
    from genomicode import jmath

    all_genes = pos_genes + neg_genes
    x = matrixlib.find_best_row_header(MATRIX, all_genes)
    header, num_found, found, not_found = x
    assert num_found, "I could not find any genes in gene set %s." % name

    MATRIX_p = MATRIX.matrix(row=pos_genes, row_header=header)
    MATRIX_n = MATRIX.matrix(row=neg_genes, row_header=header)

    num_rows = MATRIX_p.nrow() + MATRIX_n.nrow()
    if lock:
        lock.acquire()
    x = "Gene set %s contains %d genes and matched %d rows in %s."
    print x % (name, len(all_genes), num_rows, matrix_name)
    sys.stdout.flush()
    if lock:
        lock.release()

    X_p = MATRIX_p._X
    X_n = []
    for x in MATRIX_n._X:
        x = [-x for x in x]
        X_n.append(x)
    X_pn = X_p + X_n
    score = jmath.mean(X_pn, byrow=False)

    return score

def score_gene_set(gs_name, pos_genes, neg_genes, matrix_name, MATRIX,
                   lock=None):
    import arrayio
    
    # Return dict of (matrix_name, gs_name, sample) -> score.
    scores = _score_gene_set_h(
        MATRIX, matrix_name, gs_name, pos_genes, neg_genes, lock)
    sample_names = MATRIX.col_names(arrayio.COL_ID)
    assert len(sample_names) == len(scores)
    
    results = {}
    for (sample, score) in zip(sample_names, scores):
        key = matrix_name, gs_name, sample
        assert key not in results, "Duplicate: %s" % key
        results[key] = score
    return results

def score_many(jobs, lock=None):
    import arrayio

    file2matrix = {}
    
    results = {}
    for x in jobs:
        gs_name, pos_genes, neg_genes, matrix_name, matrix_file = x
        if matrix_file not in file2matrix:
            x = arrayio.read(matrix_file)
            file2matrix[matrix_file] = x
        MATRIX = file2matrix[matrix_file]
        x = score_gene_set(
            gs_name, pos_genes, neg_genes, matrix_name, MATRIX, lock=lock)
        results.update(x)
    return results

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
        "--all", dest="all_gene_sets", action="store_true", default=False,
        help="Score all gene sets in the files.")
    parser.add_argument(
        "--automatch", dest="automatch", action="store_true", default=False,
        help="Will match _UP with _DN.")
    parser.add_argument(
        "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
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
    assert args.all_gene_sets or args.gene_set, \
           "Please specify one or more gene sets to score."
    if args.num_procs < 1 or args.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")

    #if args.num_procs > 1:
    #    raise NotImplementedError, "Doesn't work.  Matrix class decorator."

    if args.libpath:
        sys.path = options.libpath + sys.path
    # Import after the library path is set.
    import time
    import multiprocessing
    from genomicode import genesetlib
    from genomicode import genepattern
    
    start_time = time.time()
    
    genepattern.fix_environ_path()

    msg = "Reading gene set file."
    if len(args.geneset_files) > 1:
        msg = "Reading gene set files."
    print msg; sys.stdout.flush()
    geneset2genes = {}  # name -> list of genes
    for filename in args.geneset_files:
        for x in genesetlib.read_genesets(filename):
            name, description, genes = x
            assert name not in geneset2genes, "Duplicate geneset: %s." % name
            geneset2genes[name] = genes

    genesets = args.gene_set
    if args.all_gene_sets:
        genesets = sorted(geneset2genes)
    if args.automatch:
        genesets = sorted(genesets)
        i = 0
        while i < len(genesets)-1:
            # Already have positive and negative.
            gs1, gs2 = genesets[i], genesets[i+1]
            ugs1, ugs2 = gs1.upper(), gs2.upper()
            if gs1.find(",") >= 0:
                i += 1
                continue
            if ugs1.endswith("_DN") and ugs2.endswith("_UP") and \
                   ugs1[:-3] == ugs2[:-3]:
                x = "%s,%s" % (gs1, gs2)
                genesets[i] = x
                del genesets[i+1]
            else:
                i += 1
    #genesets = genesets[:10]

    matrix_names = [os.path.split(x)[1] for x in expression_files]

    print "Setting up jobs."; sys.stdout.flush()
    # list of gs_name, pos_genes, neg_genes, matrix_name, matrix_file
    jobs = []
    for geneset in genesets:
        pos_gs, neg_gs = parse_geneset(geneset)
        assert pos_gs in geneset2genes, \
               "I could not find gene set: %s" % pos_gs
        if neg_gs:
            assert neg_gs in geneset2genes, \
                   "I could not find gene set: %s" % neg_gs
        gs_name = pos_gs
        if neg_gs:
            gs_name = "%s/%s" % (pos_gs, neg_gs)
            
        pos_genes = geneset2genes[pos_gs]
        neg_genes = geneset2genes.get(neg_gs, [])

        for matrix_name, matrix_file in zip(matrix_names, expression_files):
            x = gs_name, pos_genes, neg_genes, matrix_name, matrix_file
            jobs.append(x)

    # Group the jobs into batches such that jobs that use the same
    # matrix are in the same batch.
    batched_jobs = {}  # matrix_file -> list of jobs
    for i in range(len(jobs)):
        batch = jobs[i][4]
        if batch not in batched_jobs:
            batched_jobs[batch] = []
        batched_jobs[batch].append(jobs[i])
    batched_jobs = batched_jobs.values()  # list of list of jobs

    # TODO: if there are too many gene sets to score for a file, split
    # it up into multiple batches.  Don't know the tradeoff between
    # reading a file twice and calculating more gene sets.
    
    print "Scoring %d jobs." % len(jobs); sys.stdout.flush()
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    pool = multiprocessing.Pool(args.num_procs)

    scores = {}   # (matrix, geneset, sample) -> score
    results = []  # AsyncResults
    for batch in batched_jobs:
        fn_keywds = {}
        fn_keywds["lock"] = lock
        if args.num_procs == 1:
            x = score_many(batch)
            scores.update(x)
        else:
            fn_args = (batch,)
            x = pool.apply_async(score_many, fn_args, fn_keywds)
            results.append(x)
    pool.close()
    pool.join()
    for x in results:
        x = x.get()
        scores.update(x)
        
    x = [(x[0], x[2]) for x in scores]
    x = sorted({}.fromkeys(x))
    all_matrix_samples = x
    x = [x[1] for x in scores]
    x = sorted({}.fromkeys(x))
    all_genesets = x

    outhandle = open(args.outfile, 'w')
    header = ["FILE", "SAMPLE"] + all_genesets
    print >>outhandle, "\t".join(header)
    for x in all_matrix_samples:
        matrix, sample = x
        x = [scores[(matrix, x, sample)] for x in all_genesets]
        x = [matrix, sample] + x
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
    
    print "Done."
    
if __name__ == '__main__':
    main()
