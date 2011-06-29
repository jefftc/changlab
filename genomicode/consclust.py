"""

Functions:
calc_consensus
count_clusters

Classes:
ConsensusScore

"""
import os, sys

class ConsensusScore:
    """
    C          Matrix of consensus scores.
    C_flat     List of the scores on the upper diagonal.
    score2cdf  Dictionary of consensus score -> CDF.
    auc        Area under the curve.
    
    """
    def __init__(self, C, C_flat, score2cdf, auc):
        self.C = C
        self.C_flat = C_flat
        self.score2cdf = score2cdf
        self.auc = auc


def calc_consensus(X):
    """X is a gene x clustering matrix of genes that have been
    clustered multiple times.  Each column is one run of the
    clustering algorithm.  Each cell is an integer indicating the
    cluster assignment of the gene in that clustering.  It can be None
    indicating that the gene was not included in that clustering.

    Returns a ConsensusScore object.

    """
    ngene = len(X)
    nruns = len(X[0])
    assert ngene > 0, "not enough genes (%d)" % ngene
    assert nruns >= 2, "not enough clusterings (%d)" % nruns
    for x in X:
        assert len(x) == nruns

    # For each pair of genes, count the number of data sets in which
    # they occur together.
    C = _make_consensus_matrix(X)

    # Get a list of the pairwise consensus values.
    C_flat = []
    for i in range(len(C)-1):
        x = C[i][i+1:]
        C_flat.extend(x)
    TOTAL = len(C)*(len(C)-1)/2
    assert len(C_flat) == TOTAL

    # For each consensus value, calculate the cdf.
    score2cdf = {}  # dict of consensus score -> CDF
    for i, score in enumerate(sorted(C_flat, reverse=True)):
        if score in score2cdf:
            continue
        cdf = 1.0 - float(i)/len(C_flat)
        score2cdf[score] = cdf

    # Calculate the area under the curve.
    scores = sorted(score2cdf)
    n = len(scores)
    auc = 0.0
    for i in range(1, len(scores)):
        auc += (scores[i]-scores[i-1]) * score2cdf[scores[i]]

    x = ConsensusScore(C, C_flat, score2cdf, auc)
    return x

def count_clusters(X):
    """X is a gene x clustering matrix of genes that have been
    clustered multiple times.  Each column is one run of the
    clustering algorithm.  Each cell is an integer indicating the
    cluster assignment of the gene in that clustering.  It can be None
    indicating that the gene was not included in that clustering.

    Returns a tuple of (nclust, nc_mean, nc_sd) where nclust is a list
    with the same number of columns as X.  Each item is the number of
    clusters in that clustering.  nc_mean is the mean of the number of
    clusters, and nc_sd is the standard deviation.

    """
    import jmath
    
    ngene = len(X)
    nruns = len(X[0])
    if not ngene or not nruns:
        return [], 0, 0
    
    nclust = [0] * nruns
    for j in range(nruns):
        x = [X[i][j] for i in range(ngene) if X[i][j] is not None]
        nclust[j] = len({}.fromkeys(x))

    nc_mean = jmath.mean(nclust)
    nc_sd = 0
    if min(nclust) != max(nclust):
        nc_sd = jmath.stddev(nclust)
    return nclust, nc_mean, nc_sd

def _make_consensus_matrix(X):
    # Return a consensus matrix from X.  Only fills in the upper
    # diagonal of the matrix.
    # For each pair of genes, count the number of data sets in which
    # they are in the same cluster.
    ngenes, nruns = len(X), len(X[0])

    C = [[0]*ngenes for i in range(ngenes)]
    for i in range(ngenes):
        X_i = X[i]
        for j in range(i, ngenes):
            X_j = X[j]

            # Count the number of clusterings this pair of genes
            # co-occur.  May be less than the total number of
            # clusterings if different subsets of genes are used each
            # time.
            t = 0

            # Count the number of clusterings each pair of genes are
            # in the same cluster.
            c = 0
            for k in range(nruns):
                if X_i[k] is None or X_j[k] is None:
                    continue
                t += 1
                c += int(X_i[k] == X_j[k])
            if t > 0:
                C[i][j] = C[j][i] = float(c) / t
    return C


## def calc_consensus_p(dataset, NUM_SAMPLES=100):
##     import math
##     import random
##     import jmath
    
##     R = jmath.start_R()
    
##     # Calculate the consensus score for this gene set.
##     X_gset = dataset.matrix()
##     C_gset = calc_consensus_score_from_matrix(X_gset)

##     rand_nlp = [None] * NUM_SAMPLES
##     for zzz in range(NUM_SAMPLES):
##         # Randomly permute the gene set.
##         X_rand = [X_gset[i][:] for i in range(len(X_gset))]
##         for i in range(len(X_rand[0])):
##             for j in range(len(X_rand)-1, 0, -1):
##                 k = random.randint(0, j)
##                 X_rand[j][i], X_rand[k][i] = X_rand[k][i], X_rand[j][i]
##         C_rand = calc_consensus_score_from_matrix(X_rand)

##         # Calculate the p-value from the Kolmogorov-Smirnov test.
##         ow = R.options("warn")
##         R.options(warn=-1)
##         x = R.ks_test(C_gset.C_flat, C_rand.C_flat, exact=True)
##         R.options(**ow)
##         p_value = max(x["p.value"], 1E-10)
##         nlp = -math.log(p_value, 10)
##         rand_nlp[zzz] = nlp
##     mean_nlp = float(sum(rand_nlp)) / len(rand_nlp)
##     return mean_nlp
