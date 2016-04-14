#!/usr/bin/env python

import sys, os

# score_many
# score_gene_set
# score_gene
#
# match_gene_sets
# has_missing_values
# 
# _parse_geneset
# _parse_gene_names


class GeneSetScore:
    # Score for a gene set on one sample.
    def __init__(self, score, background, direction, pvalue):
        self.score = score
        self.background = background
        self.direction = direction
        self.pvalue = pvalue


class GeneScore:
    def __init__(self, score):
        self.score = score


def score_many(jobs, lock=None):
    # Return dict of (matrix_name, gs_name, index, sample) ->
    # GeneSetScore or GeneScore.
    import arrayio

    file2matrix = {}
    
    results = {}
    for x in jobs:
        (gs_name, pos_genes, neg_genes, matrix_name, matrix_file,
         any_matching_gene_sets) = x
        if matrix_file not in file2matrix:
            x = arrayio.read(matrix_file)
            file2matrix[matrix_file] = x
        MATRIX = file2matrix[matrix_file]
        assert not has_missing_values(MATRIX), \
               "Matrix %s has missing values." % matrix_name
        if pos_genes or neg_genes:
            x = score_gene_set(
                gs_name, pos_genes, neg_genes, matrix_name, MATRIX,
                any_matching_gene_sets, lock=lock)
        else:
            assert not (pos_genes == [] and neg_genes == []), \
                   "Empty gene set: %s" % gs_name
            assert pos_genes is None, "Has pos genes: %s" % gs_name
            assert neg_genes is None, "Has neg genes: %s" % gs_name
            x = score_gene(gs_name, matrix_name, MATRIX)
        # TODO: should make sure we don't overwrite previous results.
        results.update(x)
    return results


def _score_gene_set_h(MATRIX, matrix_name, name, pos_genes, neg_genes, lock):
    # Return list of GeneSetScores.
    from genomicode import genesetlib

    x = genesetlib.score_geneset(MATRIX, pos_genes, neg_genes)
    MATRIX_p, MATRIX_n, num_matches, scores, scores_bgrnd, pvalues = x

    if lock:
        lock.acquire()
    x = "Gene set %s contains %d genes (%d:%d) and matched %d rows in %s." 
    print x % (
        name, len(pos_genes)+len(neg_genes), len(pos_genes), len(neg_genes), 
        num_matches, matrix_name)
    sys.stdout.flush()
    if lock:
        lock.release()

    assert len(scores) == len(scores_bgrnd)
    directions = []
    for x in zip(scores, scores_bgrnd):
        score, background = x
        x = "UP"
        if score < background:
            x = "DOWN"
        directions.append(x)

    geneset_scores = []
    for x in zip(scores, scores_bgrnd, directions, pvalues):
        score, background, direction, pvalue = x
        x = GeneSetScore(score, background, direction, pvalue)
        geneset_scores.append(x)
    return geneset_scores


def score_gene_set(gs_name, pos_genes, neg_genes, matrix_name, MATRIX,
                   ignore_gene_not_found, lock=None):
    # Return dict of (matrix_name, gs_name, index, sample) ->
    # GeneSetScore.
    import arrayio

    try:
        scores = _score_gene_set_h(
            MATRIX, matrix_name, gs_name, pos_genes, neg_genes, lock)
    except AssertionError, x:
        if ignore_gene_not_found and \
               str(x) == "I could not find any genes in the gene set.":
            return {}
        raise
    sample_names = MATRIX.col_names(arrayio.COL_ID)
    assert len(sample_names) == len(scores)
    
    results = {}
    for i, x in enumerate(zip(sample_names, scores)):
        sample, geneset_score = x
        key = matrix_name, gs_name, i, sample
        #assert key not in results, "Duplicate: %s" % str(key)
        assert key not in results, "Duplicate: %s" % sample
        results[key] = geneset_score
    return results


def score_gene(gene_name, matrix_name, MATRIX):
    # Return dict of (matrix_name, gene_name, index, sample) ->
    # GeneScore.  Used to add gene expression to the output.
    import arrayio

    I_row, x = MATRIX._index(row=gene_name)
    if not I_row:
        return {}
    I = I_row[0]
    scores = MATRIX._X[I]
    sample_names = MATRIX.col_names(arrayio.COL_ID)
    assert len(sample_names) == len(scores)

    results = {}
    for i, (sample, score) in enumerate(zip(sample_names, scores)):
        key = matrix_name, gene_name, i, sample
        assert key not in results, "Duplicate: %s" % sample
        results[key] = GeneScore(score)
    return results


def match_gene_sets(genesets):
    # Return a list of genesets where the _UP and _DN (_DOWN) are
    # combined.
    genesets = sorted(genesets)
    i = 0
    while i < len(genesets)-1:
        gs1, gs2 = genesets[i], genesets[i+1]
        ugs1, ugs2 = gs1.upper(), gs2.upper()
        
        # If there is already a positive and negative, skip it.
        if gs1.find(",") >= 0:
            i += 1
            continue
        if ugs1.endswith("_DOWN"):
            ugs1 = ugs1[:-5] + "_DN"
        if ugs1.endswith("_DN") and ugs2.endswith("_UP") and \
               ugs1[:-3] == ugs2[:-3]:
            x = "%s,%s" % (gs2, gs1)
            genesets[i] = x
            del genesets[i+1]
        else:
            i += 1
    return genesets


def has_missing_values(MATRIX):
    for x in MATRIX._X:
        if None in x:
            return True
    return False


def _parse_geneset(geneset):
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


def _parse_gene_names(gene_name_list):
    # This can a list of comma separated genes, e.g.
    # ["E2F1", "E2F2,E2F3"]
    # Need to separate them out.
    gene_names = []
    for x in gene_name_list:
        x = x.split(",")
        gene_names.extend(x)
    return gene_names


def main():
    import argparse
    import glob
    import itertools
    
    DEF_PVALUE = 0.05

    parser = argparse.ArgumentParser(
        description="Score a gene set on a gene expression data set.")

    parser.add_argument(
        "expression_files", nargs="+", help="Data set(s) to score.")
    parser.add_argument(
        "-o", dest="outfile", default=None, help="Name of file for results.")
    parser.add_argument(
        "--transpose", action="store_true",
        help="Transpose the output matrix.")
    parser.add_argument(
        "--pvalue", type=float, default=DEF_PVALUE,
        help="p-value cutoff for determining significant changes "
        "(default %g)." % DEF_PVALUE)
    
    parser.add_argument(
        "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")

    # Assumes that there are no commas in names of gene sets.
    group = parser.add_argument_group(title="Gene Set")
    group.add_argument(
        "--geneset_file", dest="geneset_files", action="append", default=[],
        help="File(s) with gene sets.  Should be in gmx or gmt format.")
    group.add_argument(
        "-g", dest="gene_set", action="append", default=[],
        help="Name of the gene set to score.  If you want to score both "
        "the positively and negatively correlated genes, specify both "
        "gene sets using the format: <positive_geneset>,<negative_geneset>.  "
        "You can use this option multiple times to score more than one gene "
        "set.")
    group.add_argument(
        "--all", dest="all_gene_sets", action="store_true", default=False,
        help="Score all gene sets in the files.")
    group.add_argument(
        "--any_matching", dest="any_matching_gene_sets",
        action="store_true", default=False,
        help="Score gene sets in the files that matches these genes.")
    group.add_argument(
        "--automatch", action="store_true", default=False,
        help="Will match _UP with _DN (or _DOWN).")
    
    group = parser.add_argument_group(
        title="Genes", description="Add gene expression profiles to output.")
    group.add_argument(
        "--genes", default=[], action="append",
        help="Comma-separated list of IDs (e.g. probes, gene names) "
        "to include.")

    args = parser.parse_args()
    assert args.expression_files, \
           "Please specify an expression data set to score."
    expression_files = []
    for x in args.expression_files:
        xg = glob.glob(x)
        assert xg, "I could not find the expression file: %s" % x
        expression_files.extend(xg)
    for x in expression_files:
        assert os.path.exists(x), \
           "I could not find the expression file: %s" % x
    assert args.outfile, "Please specify the name of an outfile."
    
    if args.num_procs < 1 or args.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")
    assert args.pvalue > 0 and args.pvalue <= 1, \
           "Invalid pvalue %g" % args.pvalue

    assert args.geneset_files, "Please specify one or more geneset files."
    for x in args.geneset_files:
        assert os.path.exists(x), "I could not find the gene set file: %s" % x
    assert args.all_gene_sets or args.gene_set or args.any_matching_gene_sets,\
           "Please specify one or more gene sets to score."
    if args.all_gene_sets:
        assert not args.gene_set and not args.any_matching_gene_sets
    if args.any_matching_gene_sets:
        assert not args.gene_set and not args.all_gene_sets
    
    #if args.num_procs > 1:
    #    raise NotImplementedError, "Doesn't work.  Matrix class decorator."

    if args.libpath:
        sys.path = args.libpath + sys.path
    # Import after the library path is set.
    #import time
    import multiprocessing
    from genomicode import genesetlib
    from genomicode import genepattern
    from genomicode import jmath
    
    #start_time = time.time()
    
    genepattern.fix_environ_path()

    gene_names = _parse_gene_names(args.genes)
    
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
    if args.all_gene_sets or args.any_matching_gene_sets:
        genesets = sorted(geneset2genes)
    if args.automatch:
        genesets = match_gene_sets(genesets)
    #genesets = genesets[:10]

    matrix_names = [os.path.split(x)[1] for x in expression_files]

    print "Setting up jobs."; sys.stdout.flush()
    ignore_gene_not_found = args.any_matching_gene_sets
    # list of gs_name, pos_genes, neg_genes, matrix_name, matrix_file
    # list of gene_name, None, None, matrix_name, matrix_file
    jobs = []
    for geneset in genesets:
        pos_gs, neg_gs = _parse_geneset(geneset)
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

        if not pos_genes and not neg_genes:
            print "Empty gene set: %s.  Skipping." % gs_name
            continue

        for matrix_name, matrix_file in zip(matrix_names, expression_files):
            x = gs_name, pos_genes, neg_genes, matrix_name, matrix_file, \
                ignore_gene_not_found
            jobs.append(x)
    for name in gene_names:
        for matrix_name, matrix_file in zip(matrix_names, expression_files):
            x = name, None, None, matrix_name, matrix_file, None
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

    job_str = "jobs"
    if len(jobs) == 1:
        job_str = "job"
    print "Scoring %d %s." % (len(jobs), job_str); sys.stdout.flush()
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    pool = multiprocessing.Pool(args.num_procs)

    # (matrix, geneset, index, sample) -> GeneSetScore or GeneScore
    score_dict = {}
    results = []  # AsyncResults
    for batch in batched_jobs:
        fn_args = (batch,)
        fn_keywds = {}
        fn_keywds["lock"] = lock
        if args.num_procs == 1:
            x = score_many(batch)
            score_dict.update(x)
        else:
            x = pool.apply_async(score_many, fn_args, fn_keywds)
            results.append(x)
    pool.close()
    pool.join()
    for x in results:
        x = x.get()
        score_dict.update(x)

    all_matrix_samples = []
    all_genesets = []
    all_genes = []
    for (x, score) in score_dict.iteritems():
        matrix_name, gene_name, index, sample = x
        x = matrix_name, index, sample
        all_matrix_samples.append(x)
        if isinstance(score, GeneSetScore):
            all_genesets.append(gene_name)
        elif isinstance(score, GeneScore):
            all_genes.append(gene_name)
        else:
            raise AssertionError
    all_matrix_samples = sorted({}.fromkeys(all_matrix_samples))
    all_genesets = sorted({}.fromkeys(all_genesets))
    all_genes = sorted({}.fromkeys(all_genes))

    # Format the output.  Columns should be in order:
    # <SAMPLE> <FILE>
    # <GS SCORES> ... <GS DIRECTION> ... <GS PVALUE> ... <GS SIGNIFICANT> ...
    # <GENES> ...
    header = ["SAMPLE", "FILE"]
    x = ["", "direction", "pvalue", "significant"]
    for x in itertools.product(x, all_genesets):
        suffix, name = x
        x = "%s %s" % (name, suffix)
        x = x.strip()
        header = header + [x]
    for g in all_genes:
        header = header + [g]
        
    output = []
    output.append(header)
    for x in all_matrix_samples:
        matrix, index, sample = x
        #x = [scores[(matrix, x, index, sample)] for x in all_genesets]

        # Get the scores for the gene sets.
        keys = [(matrix, x, index, sample) for x in all_genesets]
        default = GeneSetScore("", "", "", "")
        scores = [score_dict.get(x, default).score for x in keys]
        directs = [score_dict.get(x, default).direction for x in keys]
        pvalues = [score_dict.get(x, default).pvalue for x in keys]
        signifs = []
        for x in zip(directs, pvalues):
            direct, pvalue = x
            x = ""
            if type(pvalue) is type(0.0) and pvalue < args.pvalue:
                x = direct
            signifs.append(x)

        # Get the scores for the genes.
        keys = [(matrix, x, index, sample) for x in all_genes]
        default = GeneScore("")
        gene_scores = [score_dict.get(x, default).score for x in keys]

        x = [sample, matrix] + \
            scores + directs + pvalues + signifs + gene_scores
        assert len(x) == len(header)
        output.append(x)

    if args.transpose:
        output = jmath.transpose(output)
        
    outhandle = open(args.outfile, 'w')
    for x in output:
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
    
    print "Done."

    
if __name__ == '__main__':
    main()
