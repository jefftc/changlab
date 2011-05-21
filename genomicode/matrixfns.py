"""

This module deprecated in favor of jmath (100315).

BUG: many functions fail when given lists.

slice    Take a slice from a matrix.
cor      Calculate the correlation coefficients between each column of X.
rank     Calculate the ranks of each column of X.
log      Take the log of a matrix.
flatten  Convert 2D matrix to list.

shapiro_test

"""

import os, sys

R = None
def start_R():
    global R
    if R is None:
        import rpy
        R = rpy.r
    return R

def cor(M, bycol=1, method=None):
    method = method or "pearson"

    if method == "pearson":
        import numpy
        # numpy calculates correlations by rows.  If want by column,
        # then transpose it.
        rowvar = 1
        if bycol:
            rowvar = 0
        return numpy.corrcoef(M, rowvar=rowvar).tolist()

    # R calculates correlations by columns.  If want by row, then
    # transpose it.
    if not bycol:
        M = numpy.transpose(M)
    R = start_R()
    X = flatten(t(M))   # flatten as column major
    X_str = ",".join(map(str, X))
    R("Y <- c(%s)" % X_str)
    R("Y <- matrix(Y, %d, %d)" % (nrow(M), ncol(M)))

    # Turn off potential warning.
    # Warning message:the standard deviation is zero in:
    # cor(x, y, na.method, method == "kendall")
    R('ow <- options("warn")')
    R('options(warn=-1)')
    cors = R('cor(Y, method="%s")' % method).tolist()
    R('options(ow)')
    return cors

    #import numpy
    #return numpy.corrcoef(numpy.transpose(M)).tolist()
    #cors = [[None]*ncol(X) for i in range(ncol(X))]
    #for i in range(ncol(X)):
    #    for j in range(i, ncol(X)):
    #        x1 = slice(X, None, i)
    #        x2 = slice(X, None, j)
    #        cors = numpy.corrcoef([x1, x2])
    #        r = cors[0][1]
    #        cors[i][j] = r
    #        cors[j][i] = r
    #return cors

def rank(M):
    # Calculate the ranks of each column.
    R = start_R()
    X = flatten(t(M))   # flatten as column major
    X_str = ",".join(map(str, X))
    R("X <- c(%s)" % X_str)
    R("X <- matrix(X, %d, %d)" % (nrow(M), ncol(M)))
    X_rank = R("apply(X, 2, rank)").tolist()
    return X_rank

def log(M, base=None):
    import math
    den = 1
    if base is not None:
        den = math.log(base)
    M_l = [[math.log(x)/den for x in row] for row in M]
    return M_l

def flatten(M):
    l = [None] * (nrow(M)*ncol(M))
    i = 0
    for x in M:
        l[i:i+len(x)] = x
        i += len(x)
    return l

def shapiro_test(data):
    # Only will work with 5000 samples.  If too many, then take a
    # random sample.
    import random

    assert len(data) >= 3, "need at least 3 samples"

    MAXLEN = 5000
    if len(data) > MAXLEN:
        data_i = {}
        while len(data_i) < MAXLEN:
            data_i[random.randint(0, len(data)-1)] = 1
        data = [data[i] for i in data_i]
        #schwartz = [(random.random(), x) for x in data]
        #schwartz.sort()
        #data = [x[-1] for x in schwartz]
        #data = data[:5000]
    R = start_R()
    x = R.shapiro_test(data)
    return x["statistic"]["W"], x["p.value"]
