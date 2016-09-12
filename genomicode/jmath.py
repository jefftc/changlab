"""

Many of the functions can be applied to single objects, lists, or
matrices.
<fn>_item
<fn>_list
<fn>_matrix


Functions:
safe_int
safe_float

is_list
is_matrix
nrow
ncol
dim

log      Take the log of a single number, a list, or a matrix.
exp
mysum
mean
safe_mean
median
safe_median
var
iqr      inter-quartile range
stddev
cov
cor
order
rank
min
max
which_max
svd
eig

log_add

equal_matrix

cmh_bonferroni   Correct multiple hypotheses, Bonferroni.
cmh_fdr          Correct multiple hypotheses, False Discovery Rate.

fisher_z
norm_mv
safe_norm_mv

transpose
flatten

dnorm
rmvnorm
dparzen
pparzen
qparzen

choose
lchoose
factorial
lfactorial
sample

int_simpson
int_simpson_adaptive

fisher_test
wilcox_test
shapiro_test

start_R
R_equals
R2py_matrix

apply
match    Return list of indexes of matches from list1 in list2.

"""
import math

class _fn:
    def __init__(self, fn, *args, **keywds):
        self.fn = fn
        self.args = args
        self.keywds = keywds
    def __call__(self, X):
        return self.fn(X, *self.args, **self.keywds)

def _dispatch(X, single_fn, list_fn, matrix_fn):
    import Matrix

    # If it's a Matrix object, convert to a list of lists.
    if isinstance(X, Matrix.AbstractMatrix):
        X = X.slice()
    
    if is_matrix(X):
        assert matrix_fn, "matrix not supported or implemented"
        return matrix_fn(X)
    elif is_list(X):
        assert list_fn, "list not supported or implemented"
        return list_fn(X)
    assert single_fn, "single object not supported or implemented"
    return single_fn(X)


def is_int(x, allow_special=False):
    if allow_special and x in ["na", "null", "nan"]:
        return True
    try:
        int(x)
    except ValueError, x:
        return False
    return True


def safe_int(x):
    if x is None:
        return None
    if type(x) is type("") and x.lower() in ["", "na", "null"]:
        return None
    if type(x) is type("") and x.lower() == "nan":
        return None
    return int(x)

def safe_float(x):
    if x is None:
        return None
    if type(x) is type("") and x.lower() in ["", "na", "null", "-"]:
        return None
    if type(x) is type("") and x.lower() == "nan":
        return float('nan')
    return float(x)

def is_list(X):
    import types
    return type(X) in [types.ListType, types.TupleType] and not is_matrix(X)

def is_matrix(X):
    import types
    return (type(X) in [types.ListType, types.TupleType] and
            X and type(X[0]) in [types.ListType, types.TupleType])

def nrow(M):
    return len(M)

def ncol(M):
    if not len(M):
        return 0
    return len(M[0])

def dim(M):
    return nrow(M), ncol(M)

def log_item(X, base=None, safe=0):
    if safe and X<=0:
        return 0.0
    den = 1
    if base is not None:
        den = math.log(base)
    return math.log(X)/den

def log_list(X, base=None, safe=0):
    #return [log_item(x, base=base, safe=safe) for x in X]
    # Optimization: write this directly into here.
    den = 1
    if base is not None:
        den = math.log(base)
    Y = [None] * len(X)
    for i in range(len(X)):
        if safe and X[i] <= 0:
            Y[i] = 0.0
        else:
            Y[i] = math.log(X[i])/den
    return Y

def log_matrix(X, base=None, safe=0):
    return [log_list(x, base=base, safe=safe) for x in X]

def log(X, base=None, safe=0):
    return _dispatch(
        X, _fn(log_item, base, safe), _fn(log_list, base, safe),
        _fn(log_matrix, base, safe))

def exp_item(X, base=None):
    # Return base raised to the power of X.
    if base is not None:
        x = math.exp(X)
    else:
        x = base ** X
    return x

def exp_list(X, base=None):
    return [exp_item(x, base=base) for x in X]

def exp_matrix(X, base=None):
    return [exp_list(x, base=base) for x in X]

def exp(X, base=None):
    return _dispatch(
        X, _fn(exp_item, base), _fn(exp_list, base), _fn(exp_matrix, base))

def sum_list(X):
    return sum(X)

def sum_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [sum(x) for x in X]

def mysum(X, byrow=1):
    return _dispatch(X, None, sum_list, _fn(sum_matrix, byrow=byrow))

def prod_list(X):
    total = 1
    for x in X:
        total = total * x
    return total

def prod(X, byrow=1):
    return _dispatch(X, None, sum_list, None)


def mean_list(X):
    assert len(X) > 0
    return float(sum(X)) / len(X)

def mean_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [mean_list(x) for x in X]

def mean(X, byrow=1):
    return _dispatch(X, None, mean_list, _fn(mean_matrix, byrow=byrow))

def safe_mean_list(X):
    import math
    assert len(X) > 0
    X = [x for x in X if x is not None and not math.isnan(x)]
    if not X:
        return 0
    return float(sum(X)) / len(X)

def safe_mean_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [safe_mean_list(x) for x in X]

def safe_mean(X, byrow=1):
    return _dispatch(
        X, None, safe_mean_list, _fn(safe_mean_matrix, byrow=byrow))

def median_list(X):
    assert len(X) > 0
    X = sorted(X)
    if len(X) % 2:
        # odd number of items.
        i = len(X)/2
        return X[i]
    i2 = len(X)/2
    i1 = i2-1
    return (X[i1]+X[i2])/2.0

def median_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [median_list(x) for x in X]
    
def median(X, byrow=1):
    return _dispatch(X, None, _fn(median_list), _fn(median_matrix,byrow=byrow))

def safe_median_list(X):
    import math
    assert len(X) > 0
    X = [x for x in X if x is not None and not math.isnan(x)]
    if not X:
        return 0
    return median_list(X)

def safe_median_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [safe_median_list(x) for x in X]

def safe_median(X, byrow=1):
    return _dispatch(
        X, None, _fn(safe_median_list), _fn(safe_median_matrix, byrow=byrow))

def var_list(X):
    assert len(X) > 0
    m = mean_list(X)
    total = sum([(x-m)**2 for x in X])
    return float(total) / (len(X)-1)

def var_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [var_list(x) for x in X]

def var(X, byrow=1):
    return _dispatch(X, None, _fn(var_list), _fn(var_matrix, byrow=byrow))

def safe_var_list(X):
    assert len(X) > 0
    X = [x for x in X if x is not None]
    if len(X) <= 1:
        return 0
    return var_list(X)

def safe_var_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [safe_var_list(x) for x in X]

def safe_var(X, byrow=1):
    return _dispatch(
        X, None, _fn(safe_var_list), _fn(safe_var_matrix, byrow=byrow))

def iqr_list(X):
    assert len(X) > 0

    # 75% - 25%
    X = sorted(X)
    i25 = int(round((len(X)-1) * 0.25))
    i75 = int(round((len(X)-1) * 0.75))
    assert i25 >= 0 and i25 < len(X)
    assert i75 >= 0 and i75 < len(X)
    q25 = X[i25]
    q75 = X[i75]
    return q75 - q25

def iqr_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [iqr_list(x) for x in X]

def iqr(X, byrow=1):
    return _dispatch(X, None, _fn(iqr_list), _fn(iqr_matrix, byrow=byrow))

def stddev_list(X):
    v = var(X)
    return math.sqrt(v)

def stddev_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    return [stddev_list(x) for x in X]

def stddev(X, byrow=1):
    return _dispatch(X, None, _fn(stddev_list),_fn(stddev_matrix, byrow))

def cov_list(X, Y):
    assert Y is not None
    M = [X, Y]
    x = cov_matrix(M, byrow=1)
    return x[0][1]

def cov_matrix(X, byrow=1):
    import numpy
    return numpy.cov(X, None, rowvar=byrow, bias=0).tolist()

def cov(X, Y=None, byrow=1):
    return _dispatch(X, None, _fn(cov_list, Y), _fn(cov_matrix, byrow=byrow))

def cor_list(X, Y, safe=0, abs=0):
    assert Y is not None
    M = [X, Y]
    x = cor_matrix(M, byrow=1, safe=safe, abs=abs)
    return x[0][1]

def cor_matrix(X, byrow=1, safe=0, abs=0):
    import numpy
    
    if not byrow:
        raise NotImplementedError

    if safe and ncol(X) == 1:
        # If there's only 1 sample, there will be a divide-by-zero in
        # cov.
        cor = [[0]*nrow(X) for i in range(nrow(X))]
        return cor
    cor = cov(X, byrow=byrow)

    sd = numpy.sqrt(numpy.diag(cor))
    for i in range(nrow(cor)):
        for j in range(i, ncol(cor)):
            if safe and (not sd[i] or not sd[j]):
                # Ignore (zero-out) the genes with no variance.
                cor[i][j] = 0.0
                cor[j][i] = 0.0
            else:
                cor[i][j] = cor[i][j] / (sd[i]*sd[j])
                cor[j][i] = cor[i][j]
    return cor

## def cor(M, bycol=1, method=None):
##     method = method or "pearson"

##     if method == "pearson":
##         import numpy
##         # numpy calculates correlations by rows.  If want by column,
##         # then transpose it.
##         rowvar = 1
##         if bycol:
##             rowvar = 0
##         return numpy.corrcoef(M, rowvar=rowvar).tolist()

##     # R calculates correlations by columns.  If want by row, then
##     # transpose it.
##     if not bycol:
##         M = numpy.transpose(M)
##     R = start_R()
##     X = flatten(t(M))   # flatten as column major
##     X_str = ",".join(map(str, X))
##     R("Y <- c(%s)" % X_str)
##     R("Y <- matrix(Y, %d, %d)" % (nrow(M), ncol(M)))

##     # Turn off potential warning.
##     # Warning message:the standard deviation is zero in:
##     # cor(x, y, na.method, method == "kendall")
##     R('ow <- options("warn")')
##     R('options(warn=-1)')
##     cors = R('cor(Y, method="%s")' % method).tolist()
##     R('options(ow)')
##     return cors

##     #import numpy
##     #return numpy.corrcoef(numpy.transpose(M)).tolist()
##     #cors = [[None]*ncol(X) for i in range(ncol(X))]
##     #for i in range(ncol(X)):
##     #    for j in range(i, ncol(X)):
##     #        x1 = slice(X, None, i)
##     #        x2 = slice(X, None, j)
##     #        cors = numpy.corrcoef([x1, x2])
##     #        r = cors[0][1]
##     #        cors[i][j] = r
##     #        cors[j][i] = r
##     #return cors

def cor(X, Y=None, byrow=1, safe=0, abs=0):
    return _dispatch(
        X, None, _fn(cor_list, Y, safe=safe, abs=abs),
        _fn(cor_matrix, byrow=byrow, safe=safe, abs=abs))

def order_list(X, decreasing=0):
    I = range(len(X))
    schwartz = [(x, i) for (x, i) in zip(X, I)]
    schwartz.sort()
    I = [x[-1] for x in schwartz]
    if decreasing:
        I.reverse()
    return I

def order_matrix(X, byrow=1, decreasing=0):
    if not byrow:
        X = transpose(X)
    X_order = [order_list(x) for x in X]
    if not byrow:
        X_order = transpose(X_order)
    return X_order

def order(X, byrow=1, decreasing=0):
    return _dispatch(
        X, None, _fn(order_list, decreasing=decreasing),
        _fn(order_matrix, byrow=byrow, decreasing=decreasing))

def rank_list(X):
    # Return the ranks, where the smallest value of X is ranked 0, and
    # the largest is ranked n-1.
    value2count = {}
    for x in X:
        if x not in value2count:
            value2count[x] = 1
        else:
            value2count[x] += 1
    value2rank = {}
    rank = 0
    for x in sorted(value2count):
        count = value2count[x]
        if count == 1:
            value2rank[x] = float(rank)
            rank += 1
        else:
            value2rank[x] = rank + (count-1)/2.0
            rank += count
    X_rank = [value2rank[x] for x in X]
    return X_rank
    #R = start_R()
    #X_str = ",".join(map(str, X))
    #R("x <- c(%s)" % X_str)
    #X_rank = R("rank(x)")
    ## R is 1-based, and we want 0-based rankings.
    #X_rank = [x-1 for x in X_rank]
    #return X_rank

def rank_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    X_rank = [rank_list(x) for x in X]
    if not byrow:
        X_rank = transpose(X_rank)
    return X_rank

def rank(X, byrow=1):
    return _dispatch(X, None, _fn(rank_list), _fn(rank_matrix, byrow=byrow))

max_ = max
min_ = min

def min_list(X):
    return min_(X)

def min_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    X_min = [min_list(x) for x in X]
    return X_min

def min(X, byrow=1):
    return _dispatch(X, None, _fn(min_list), _fn(min_matrix, byrow=byrow))

def max_list(X):
    return max_(X)

def max_matrix(X, byrow=1):
    if not byrow:
        X = transpose(X)
    X_max = [max_list(x) for x in X]
    return X_max

def max(X, byrow=1):
    return _dispatch(X, None, _fn(max_list), _fn(max_matrix, byrow=byrow))

def which_max(X):
    max_val = max_i = None
    for i in range(len(X)):
        if max_val is None or X[i] > max_val:
            max_val, max_i = X[i], i
    return max_i

def svd(X):
    # Return U, S, Vt.
    # U   nrow x k  columns are principal components
    # S   k array
    # Vt  ncol x k  columns are principal components
    # k   min(nrow, ncol)
    import numpy
    
    # U  nrow x k  columns are principal components
    # V  k x ncol  rows are principal components
    U, s, V = numpy.linalg.svd(X, full_matrices=False)
    U = U.tolist()
    s = s.tolist()
    V = numpy.transpose(V).tolist()
    return U, s, V


def eig(X):
    # Return:
    # list of eigenvalues, sorted in decreasing order
    # matrix of right eigenvectors (each column is a vector)
    # X x = lambda x
    import numpy

    # Bug: should check if X is square.
    w, v = numpy.linalg.eig(X)
    w = w.tolist()
    v = v.tolist()

    # For complex numbers, just sort by the real part.
    w_real = w
    if type(w[0]) is type(complex()):
        w_real = [x.real for x in w]
    
    O = order(w_real, decreasing=1)
    vt = transpose(v)   # transpose from columns to rows
    w = [w[i] for i in O]
    vt = [vt[i] for i in O]
    v = transpose(vt)
    
    return w, v
    

def log_add(logx, logy):
    # return log(x+y), given log(x) and log(y), checking for overflow
    # problems
    import math
    if logy - logx > 100:
        return logy
    elif logx - logy > 100:
        return logx
    minxy = min_(logx, logy)
    return minxy + math.log(math.exp(logx-minxy) + math.exp(logy-minxy))

def equal_matrix(X1, X2, precision=1E-8):
    if not X1 and not X2:
        return True
    assert len(X1) == len(X2)
    assert len(X1[0]) == len(X2[0])

    for i in range(len(X1)):
        for j in range(len(X1[i])):
            if abs(X1[i][j]-X2[i][j]) > precision:
                return False
    return True

def assert_equal_matrix(X1, X2, precision=1E-6):
    if not X1 and not X2:
        return
    assert len(X1) == len(X2)
    assert len(X1[0]) == len(X2[0])

    max_diff = 0
    for i in range(len(X1)):
        for j in range(len(X1[i])):
            diff = abs(X1[i][j]-X2[i][j])
            max_diff = max_(diff, max_diff)
            if diff > precision:
                raise AssertionError, "%d, %d: %s %s" % (
                    i, j, X1[i][j], X2[i][j])
    return max_diff

def cmh_bonferroni(p_values):
    # p_values should be a list of float or None.  Will return a
    # parallel list of Bonferroni corrected p-values or None.

    # Check types of p_values.
    for x in p_values:
        assert type(x) in [type(None), type(0.0), type(0)]

    x = [x for x in p_values if x is not None]
    m = len(x)
    bonf = [None] * len(p_values)
    for i in range(len(p_values)):
        if p_values[i] is None:
            continue
        bonf[i] = min_(p_values[i]*m, 1)
    #x = [p*m for p in p_values]
    #bonf = [min(x, 1) for x in x]
    return bonf
    
def cmh_fdr_bh(p_values):
    # Check types of p_values.
    for x in p_values:
        assert type(x) in [type(None), type(0.0), type(0)]
    p_clean = []
    p_I = []
    for i, p in enumerate(p_values):
        if p is None:
            continue
        p_I.append(i)
        p_clean.append(p)
    
    m = len(p_clean)

    # Sort p.values in increasing order.
    O = order(p_clean)
    I = range(1, m+1)
    O_rev = order([I[o] for o in O])
    p_clean = [p_clean[o] for o in O]

    k = range(1, m+1)
    x = [float(m)/k[i] * p_clean[i] for i in range(len(p_clean))]
    x = [min_(x, 1) for x in x]
    x = [min_(x[i:]) for i in range(len(x))]
    fdr_clean = [x[o] for o in O_rev]


    fdr = [None] * len(p_values)
    for f, i in zip(fdr_clean, p_I):
        fdr[i] = f
    return fdr

def fisher_z_item(R, N):
    # Fisher transformation.  N is sample size.
    # Need at least N=2, or will generate divide-by-zero error.
    assert N >= 2
    assert R is not None
    assert R >= -1 and R <= 1

    EPS = 1E-10
    R = min_(R, 1-EPS)
    R = max_(R, -1+EPS)

    if N <= 25:
        # For Hotelling's transform.
        var_W = 1.0 / (N-1)
    else:
        var_W = 1.0 / (N-3)

    W = 0.5 * math.log((1+R)/(1-R))
    if N <= 25:
        # Use Hotelling's transform.
        W = W - (3*W+math.tanh(W)) / (4*(N-1))
    z = W / math.sqrt(var_W)
    return z

def fisher_z_list(R, N):
    return [fisher_z_item(x, N) for x in R]

def fisher_z_matrix(cor, N):
    # cor is a matrix of correlations.
    # Do Fisher's transform on correlations to yield z-scores.
    zcors = [[0]*ncol(cor) for i in range(nrow(cor))]
    MAX, MIN = 1-1E-10, -1+1E-10
    for i in range(nrow(cor)):
        for j in range(i, ncol(cor)):
            r = cor[i][j]
            #r0 = max(min(r, 1-1E-10), -1+1E-10)    # slow
            r0 = r
            if r0 > MAX:
                r0 = MAX
            elif r0 < MIN:
                r0 = MIN

            z = _fisher_z_1d(r0, N)
            zcors[i][j] = z
            zcors[j][i] = z
    return zcors

def fisher_z(R, N):
    # Apply a Fisher's Z transform to either a correlation or a matrix
    # of correlations.  N is the number of samples that were
    # correlated.
    return _dispatch(
        R, _fn(fisher_z_item, N), _fn(fisher_z_list, N),
        _fn(fisher_z_matrix, N))

def norm_mv_list(x, M=0, V=1):
    V_0 = var(x)
    M_0 = mean(x)
    x = [(x-M_0)*math.sqrt(V/V_0)+M for x in x]
    return x

def norm_mv_matrix(X, M=0, V=1, byrow=1):
    if not byrow:
        X = transpose(X)
    X_norm = [norm_mv_list(x, M=M, V=V) for x in X]
    if not byrow:
        X_norm = transpose(X_norm)
    return X_norm

def norm_mv(X, M=0, V=1, byrow=1):
    # Normalize a list of numbers so the mean is M and variance is V.
    return _dispatch(
        X, None, _fn(norm_mv_list, M, V), _fn(norm_mv_matrix, M, V, byrow))

def safe_norm_mv_list(X, M=0, V=1):
    assert len(X) > 0
    V_0 = safe_var(X)
    M_0 = safe_mean(X)
    X_norm = X[:]
    for i in range(len(X)):
        if X[i] is None:
            continue
        if abs(V_0) < 1E-100:
            X_norm[i] = M
        else:
            X_norm[i] = (X[i]-M_0)*math.sqrt(V/V_0)+M
    return X_norm

def safe_norm_mv_matrix(X, M=0, V=1, byrow=1):
    if not byrow:
        X = transpose(X)
    X_norm = [safe_norm_mv_list(x, M=M, V=V) for x in X]
    if not byrow:
        X_norm = transpose(X_norm)
    return X_norm

def safe_norm_mv(X, M=0, V=1, byrow=1):
    # Normalize a list of numbers so the mean is M and variance is V.
    return _dispatch(
        X, None, _fn(safe_norm_mv_list, M, V),
        _fn(safe_norm_mv_matrix, M, V, byrow))

def transpose(X):
    # Bug: Need to check dimensions of X.
    # Bug: Need to check type of X.
    nrow, ncol = len(X), len(X[0])
    t_X = [[None] * nrow for i in range(ncol)]
    for i in range(nrow):
        for j in range(ncol):
            t_X[j][i] = X[i][j]
    return t_X

def flatten(X):
    X_flat = []
    for x in X:
        if type(x) in [type(()), type([])]:
            X_flat.extend(flatten(x))
        else:
            X_flat.append(x)
    return X_flat

## def flatten(M):
##     l = [None] * (nrow(M)*ncol(M))
##     i = 0
##     for x in M:
##         l[i:i+len(x)] = x
##         i += len(x)
##     return l

SQRT_TWO_PI = math.sqrt(2 * math.pi)
def dnorm_item(X, mean=0, sd=1):
    assert sd >= 0
    a = sd * SQRT_TWO_PI
    b = X-mean
    c = -(b**2)/(2.0*sd**2)
    return 1.0/a*math.exp(c)

def dnorm_list(X, mean=0, sd=1):
    return [dnorm_item(x, mean=mean, sd=sd) for x in X]

def dnorm_matrix(X, mean=0, sd=1):
    return [dnorm_list(x, mean=mean, sd=sd) for x in X]
    
def dnorm(X, mean=0, sd=1):
    return _dispatch(
        X, _fn(dnorm_item, mean, sd), _fn(dnorm_list, mean, sd),
        _fn(dnorm_matrix, mean, sd))

def rmvnorm(n, mean=None, sigma=None):
    """Sample from a multivariate normal distribution.  mean

    """
    assert not (mean is None and sigma is None)
    if mean is None:
        mean = [0] * len(sigma)
    if sigma is None:
        sigma = [[0]*len(mean) for i in range(len(mean))]
        for i in range(len(mean)):
            sigma[i][i] = 1

    # Should check if mean is a list.
    assert is_matrix(sigma), "sigma must be a matrix"
    assert nrow(sigma) == ncol(sigma), "sigma must be a square matrix"
    assert nrow(sigma) == len(mean), "sigma and mean are unaligned"
            
    R = start_R()
    R('library("mvtnorm")')
    R_equals(mean, "m")
    R_equals(sigma, "s")
    x = R('rmvnorm(%d, m, s)' % n)
    x = list(x)
    return x

def dparzen(x, X_obs, X_count, h=1):
    # Estimate the density at x.
    # X.obs is the observed values.
    # X.count are the number of times each one is observed.
    # Uses a guassian kernel.
    N = sum(X_count)
    h = float(h)
    total = 0.0
    for x_i, x_c in zip(X_obs, X_count):
        y = (x-x_i)/h
        total += x_c * math.exp(-(y*y)/2.0)
    total = total / SQRT_TWO_PI
    return total / (N*h)

def pparzen(x, X_obs, X_count, h=1, eps=1E-10):
    # Given x, what is the density <= x.
    # Calculate a lower bound within error of eps.
    # x_i = X_obs[i]
    # c_i = X_count[i]
    # n   = len(X_obs)
    # N   = sum(X_count)
    # c_max = max(X_count)
    # x_min = min(X_obs)
    # EPS = 1/(N*h) SUM_i^n [ 1/(sqrt(2*PI)) exp(-1/2 ((x-x_i)/h)**2)  * c_i ]
    # EPS*N*h*sqrt(2*PI) = SUM_i^n [ c_i * exp(-1/2 ((x-x_i)/h)**2) ]
    # A = EPS*N*h*sqrt(2*PI)/c_max
    # A = exp(-1/2 ((x-x_min)/h)**2)
    # -2*log(A) = ((x-x_min)/h)**2
    # -sqrt(-2*log(A)) = (x-x_min)/h     # Want the lower bound
    # x_min-h*sqrt(-2*log(A)) = x
    #print "PPARZEN"
    N = sum(X_count)
    A = eps*N*h*math.sqrt(2*math.pi)/max_(X_count)
    B = h*math.sqrt(-2*math.log(A))
    a = min_(X_obs) - B
    b = x
    #print a, b, h, eps
    #print a, b, X_obs, X_count
    return _pparzen_ul(a, b, X_obs, X_count, h=h, eps=eps)

def _pparzen_ul(p_l, p_u, X_obs, X_count, h=1, eps=1E-10):
    # Integrate using composite Simpson's rule.
    f = lambda x: dparzen(x, X_obs, X_count, h)
    a, b = p_l, p_u
    c = (a+b)/2.0
    f_a = f(a); f_b = f(b); f_c = f(c)
    s = (b-a)/6.0*(f_a + 4*f_c + f_b)
    
    sum = 0.0
    stack = [(a, c, b, f_a, f_c, f_b, eps, s)]
    while stack:
        # a   b   c   d   e
        a, c, e, f_a, f_c, f_e, eps, whole = stack.pop()
        b = (a+c)/2.0; d = (c+e)/2.0
        f_b = f(b); f_d = f(d)
        left = (c-a)/6.0 * (f_a + 4*f_b + f_c)
        right = (e-c)/6.0 * (f_c + 4*f_d + f_e)
        if abs(left + right - whole) <= 15*eps:
            sum += left+right+(left+right-whole)/15
            continue
        stack.append((a, b, c, f_a, f_b, f_c, eps/2.0, left))
        stack.append((c, d, e, f_c, f_d, f_e, eps/2.0, right))
    return sum

def qparzen(q, X_obs, X_count, h=1, eps=1E-10):
    # Given a quantile, what is the value x.
    # Could optimize: if already calculated the distribution in a
    # previous call, don't recalculate that portion of the
    # distribution.  pparzen will need to take a lower and upper
    # bound.
    N = sum(X_count)
    A = eps*N*h*math.sqrt(2*math.pi)/max_(X_count)
    B = h*math.sqrt(-2*math.log(A))
    # Search over the range of all the observed values.
    x_min, x_max = min_(X_obs)-B, max_(X_obs)+B

    q_cur = 0.0
    while abs(q_cur-q) > eps:
        x = (x_min+x_max)/2.0
        d = _pparzen_ul(x_min, x, X_obs, X_count, h=h, eps=eps)
        if d <= eps:
            break
        if q_cur + d < q:
            q_cur += d
            x_min = x
        else:
            x_max = x
    return x

## def qparzen(q, X_obs, X_count, h=1, eps=1E-10):
##     # Given a quantile, what is the value x.
##     # Could optimize: if already calculated the distribution in a
##     # previous call, don't recalculate that portion of the
##     # distribution.  pparzen will need to take a lower and upper
##     # bound.
##     N = sum(X_count)
##     A = eps*N*h*math.sqrt(2*math.pi)/max(X_count)
##     B = h*math.sqrt(-2*math.log(A))
##     # Search over the range of all the observed values.
##     x_min, x_max = min(X_obs)-B, max(X_obs)+B
##     x_delta = 1E-10
##     while x_min < x_max:
##         x = (x_min+x_max)/2.0
##         q_est = pparzen(x, X_obs, X_count, h=h, eps=eps)
##         if abs(q_est - q) < eps:
##             break
##         if q_est < q:
##             # Estimate quantile too low, so raise x.
##             x_min = x + x_delta
##         else:
##             x_max = x - x_delta
##     return x

def choose(n, k):
    import gmpy
    return int(gmpy.comb(n, k))
    #R = start_R()
    #x = R("choose(%d, %d)" % (n, k))
    #x = list(x)
    #return x[0]

def lchoose(n, k):
    import math
    return math.log(choose(n, k))
    #R = start_R()
    #x = R("lchoose(%d, %d)" % (n, k))
    #x = list(x)
    #return x[0]

def factorial(x):
    import gmpy
    return int(gmpy.fac(x))
    #R = start_R()
    #x = R("factorial(%d)" % x)
    #x = list(x)
    #return x[0]

def lfactorial(x):
    R = start_R()
    x = R("lfactorial(%d)" % x)
    x = list(x)
    return x[0]

def sample(X, n, replace=False, prob=None):
    # Choose n items from X.  X must be a list of numbers.
    if not replace:
        assert n <= len(X)
    R = start_R()
    R_equals(X, "X")
    R_equals(n, "size")
    R('replace <- FALSE')
    if replace:
        R('replace <- TRUE')
    R('prob <- NULL')
    if prob is not None:
        R_equals(prob, "prob")
    x = R('sample(X, size, replace=replace, prob=prob)')
    x = list(x)
    # Always return a list.
    #if n == 1:
    #    x = [x]
    return x

def int_simpson(f, a, b):
    c = (a+b)/2.0
    h3 = (b-a)/6.0
    s = h3*(f(a) + 4*f(c) + f(b))
    return s

def int_simpson_adaptive(f, a, b, eps):
    # Optimized to minimize function calls.
    c = (a+b)/2.0
    f_a = f(a); f_b = f(b); f_c = f(c)
    s = (b-a)/6.0*(f_a + 4*f_c + f_b)
    
    sum = 0.0
    stack = [(a, c, b, f_a, f_c, f_b, eps, s)]
    while stack:
        # a   b   c   d   e
        a, c, e, f_a, f_c, f_e, eps, whole = stack.pop(0)
        b = (a+c)/2.0; d = (c+e)/2.0
        f_b = f(b); f_d = f(d)
        left = (c-a)/6.0 * (f_a + 4*f_b + f_c)
        right = (e-c)/6.0 * (f_c + 4*f_d + f_e)
        if abs(left + right - whole) <= 15*eps:
            sum += left+right+(left+right-whole)/15
            continue
        stack.append((a, b, c, f_a, f_b, f_c, eps/2.0, left))
        stack.append((c, d, e, f_c, f_d, f_e, eps/2.0, right))
    return sum

def fisher_test(x00, x01, x10, x11):
    # Return a p-value.
    R = start_R()
    #x = array.array("i", [x00, x10, x01, x11])
    x = robjects.FloatVector([x00, x10, x01, x11])
    x = R["fisher.test"](R.matrix(x, 2, 2))
    return x.rx2("p.value")[0]

def wilcox_test(X, Y, exact=True):
    R = start_R()
    x = R.wilcox_test(X, Y, exact=exact)
    return x["statistic"]["W"], x["p.value"]

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
    #x = R.shapiro_test(data)
    #return x["statistic"]["W"], x["p.value"]
    data_r = robjects.FloatVector(data)
    x = R["shapiro.test"](data_r)
    p_value = list(x.rx2("p.value"))[0]
    W = list(x.rx2("statistic").rx2("W"))[0]
    return W, p_value

R = robjects = None
def start_R():
    global R, robjects
    import os

    # No.  This doesn't seem to have any effect.
    # Set the LD_LIBRARY_PATH.  Somehow the R library path doesn't get
    # set up correctly when using rpy2.
    #PATHS = [
    #    "/usr/local/lib64/R/lib",
    #    "/usr/local/lib64",
    #    "/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.34.x86_64/jre/lib/amd64/server"
    #    ]
    #
    #LD_LIBRARY_PATH = os.environ.get("LD_LIBRARY_PATH", "")
    #paths = LD_LIBRARY_PATH.split(":")
    #for x in PATHS:
    #    if x not in paths:
    #        paths.append(x)
    #paths = [x for x in paths if x and os.path.exists(x)]
    #LD_LIBRARY_PATH = ":".join(paths)
    #os.environ["LD_LIBRARY_PATH"] = LD_LIBRARY_PATH
    
    if R is None:
        #import rpy
        #R = rpy.r
        import rpy2.robjects as robjects
        robjects = robjects   # So pychecker won't complain.
        R = robjects.r
    return R

def _fmt_R_var(var):
    if type(var) is type(""):
        if var.startswith("RVAR:"):
            var = var[5:]
        else:
            # Escape the quotes.
            var = var.replace('"', '\\"')
            var = '"%s"' % var
    elif type(var) in [type([]), type(())]:
        x = [_fmt_R_var(x) for x in var]
        var = "c(%s)" % ", ".join(x)
    var = str(var)
    return var

def R_equals_item(x, varname):
    x = _fmt_R_var(x)
    R("%s <- %s" % (varname, x))
    return varname

def R_equals_vector(X, varname):
    # Should use FloatVector instead of string converstion.
    X = [_fmt_R_var(x) for x in X]
    X_str = ",".join(map(str, X))
    R("%s <- c(%s)" % (varname, X_str))
    return varname

def R_equals_matrix(M, varname, by_row=True):
    # Bug: Only works with floats.
    temp_varname = "JMATH.R.TMP"
    R = start_R()
    x = flatten(M)
    # Figure out which values are missing, and substitute them with
    # 0.0 so FloatVector doesn't complain.
    I = [i for i in range(len(x)) if x[i] == None]
    for i in I:
        x[i] = 0.0
    x = robjects.FloatVector(x)
    # Fill in the missing values.
    for i in I:
        x[i] = robjects.NA_Real
    x = robjects.r["matrix"](x, nrow=nrow(M), ncol=ncol(M), byrow=by_row)

    # Set to a temporary variable first.  varname may not be a
    # variable, e.g. "OBJ$member".
    robjects.globalenv[temp_varname] = x
    R("%s <- %s" % (varname, temp_varname))
    #X = flatten(M)
    #X_str = ",".join(map(str, X))
    #R("%s <- c(%s)" % (varname, X_str))
    #R("%s <- matrix(%s, %d, %d)" % (varname, varname, nrow(M), ncol(M)))
    return varname

def R_equals(X, varname):
    # Should rename to R_setvalue or something.
    # Also should re-order the variables.
    return _dispatch(
        X, _fn(R_equals_item, varname), _fn(R_equals_vector, varname),
        _fn(R_equals_matrix, varname))

def R_var(name):
    return "RVAR:%s" % name

def R_fn(fn_name, *args, **keywds):
    retval = keywds.get("RETVAL")
    if "RETVAL" in keywds:
        del keywds["RETVAL"]
    R = start_R()
    params = [_fmt_R_var(x) for x in args]
    for key, value in keywds.iteritems():
        value = _fmt_R_var(value)
        params.append("%s=%s" % (key, value))
    cmd = "%s(%s)" % (fn_name, ", ".join(params))
    if retval:
        cmd = "%s <- %s" % (retval, cmd)
    #print cmd
    R(cmd)

def R2py_matrix(m):
    # m should be a rpy2.robjects.Matrix object.
    import time
    start = time.time()
    # Fastest implementation (1.2s for 37,632x2 matrix)
    pym = [[None]*m.ncol for i in range(m.nrow)]
    z = 0
    #for i in range(m.nrow):
    #    for j in range(m.ncol):
    #        pym[i][j] = m[z]
    #        z += 1
    for j in range(m.ncol):
        for i in range(m.nrow):
            pym[i][j] = m[z]
            z += 1
    # This implementation is slow (5.0s).
    #pym = []
    #for i in range(m.nrow):
    #    row = m.rx(i+1, True)
    #    x = [row[j] for j in range(m.ncol)]
    #    pym.append(x)
    # This implementation is slowest (8.5s).
    #pym = [[None]*m.ncol for i in range(m.nrow)]
    #for i in range(m.nrow):
    #    for j in range(m.ncol):
    #        pym[i][j] = m.rx(i+1, j+1)[0]
    end = time.time()
    #print end - start
    return pym

def apply(M, margin, fn):
    assert margin in [1, 2]
    if margin == 1:
        results = [fn(slice(M, i, None)) for i in range(nrow(M))]
    else:
        results = [fn(slice(M, None, i)) for i in range(ncol(M))]
    return results

def match(list1, list2):
    # Return list of indexes of matches from list1 in list2.  An index
    # will be None if that item is missing from list2.
    x2i = {}
    for i, x in enumerate(list2):
        if x not in x2i:   # save the first match only.
            x2i[x] = i
    # Not faster, more confusing.
    #for i in range(len(list2)-1, -1, -1):
    #    x = list2[i]
    #    x2i[x] = i
    return [x2i.get(x) for x in list1]

def slice(M, rows, cols):
    """Take a slice of the matrix M.  rows and cols can be:
    - single numbers
    - tuples of (start, end)
    - lists of indexes

    Single numbers and tuples have the same semantics as Python
    slices.  Lists of indexes must refer to valid indexes (including
    negative indexes).

    If both rows and cols are single numbers, then returns a single
    cell.  If exactly one of rows and cols are tuples or lists, then
    returns a list.  Otherwise, returns a matrix.

    Negative slices are allowed.

    rows can be None, which is shorthand for (0, nrow(M)).
    rows can be (start, None), which is shorthand for (start, nrow(M)).
    These also apply for cols.

    """
    import operator

    def _normalize_index(i, max_i):
        assert type(i) is type(0)
        
        # Convert negative index into positive one.
        if i < 0:
            i += max_i
        assert i >= 0 and i < max_i, "matrix index out of range"
        return [i]

    def _normalize_list(li, max_i):
        return [_normalize_index(i, max_i)[0] for i in li]
    
    def _normalize_one_slice(i, max_i):
        assert type(i) is type(0)

        if i < 0:
            i += max_i
        if i > max_i:
            i = max_i
        return i

    def _normalize_slice(i0, i1, max_i):
        if i0 is None:
            i0 = 0
        if i1 is None:
            i1 = max_i
        i0 = _normalize_one_slice(i0, max_i)
        i1 = _normalize_one_slice(i1, max_i)

        if i0 > i1:
            i0 = i1 = 0
        return range(i0, i1)

    row_is_slice = rows is None or operator.isSequenceType(rows)
    col_is_slice = cols is None or operator.isSequenceType(cols)

    if rows is None:
        rows_list = range(nrow(M))
    elif type(rows) is type(()):
        assert len(rows) == 2
        rows_list = _normalize_slice(rows[0], rows[1], nrow(M))
    elif type(rows) is type([]):
        rows_list = _normalize_list(rows, nrow(M))
    else:
        rows_list = _normalize_index(rows, nrow(M))

    if cols is None:
        cols_list = range(ncol(M))
    elif type(cols) is type(()):
        assert len(cols) == 2
        cols_list = _normalize_slice(cols[0], cols[1], ncol(M))
    elif type(cols) is type([]):
        cols_list = _normalize_list(cols, ncol(M))
    else:
        cols_list = _normalize_index(cols, ncol(M))

    # First, make the slice as a matrix.  Fix the type later.
    S = [M[i] for i in rows_list]
    S = [[x[i] for i in cols_list] for x in S]
    if not rows_list or not cols_list:
        S = [[]]

    # If both are slices, then don't change anything.
    if row_is_slice and col_is_slice:
        pass
    # If neither are slices, then convert to a single number.
    elif not row_is_slice and not col_is_slice:
        S = S[0][0]
    # If the row is a slice, (and col is single number), then return a
    # list of the column data.
    elif row_is_slice:
        S = [x[0] for x in S]
    else:
        S = S[0]
    return S

def test_slice():
    M = [[0, 1, 2], [3, 4, 5]]

    tests = [
        (M, 0, 0),                  # 0
        (M, 1, 2),                  # 5
        (M, -1, 0),                 # 3
        (M, -1, -2),                # 4
        (M, -5, 0),                 # exception
        (M, 1, 3),                  # exception
        (M, 0, (0, -1)),            # [0, 1]
        (M, 0, (-1, 0)),            # []
        (M, (-1, 0), (0, 1)),       # [[]]
        (M, (-1, 0), (-1, 0)),      # [[]]
        (M, (0, 2), (-1, 0)),       # [[]]
        (M, (0, 2), 0),             # [0, 3]
        (M, None, 0),               # [0, 3]
        (M, None, None),            # [[0, 1, 2], [3, 4, 5]]
        (M, (None, None), None),    # [[0, 1, 2], [3, 4, 5]]
        (M, (None, 1), None),       # [[0, 1, 2]]
        (M, 1, [0, 2]),             # [3, 5]
        (M, [1, 0], 0),             # [3, 0]
        (M, [1, 0], [2, 0]),        # [[5, 3], [2, 0]]
        ]

    #tests = [
    #    ]
    
    for x in tests:
        try:
            x = slice(*x)
        except Exception, x:
            print "exception"
            #raise
        else:
            print x


    import numpy
    return numpy.transpose(M).tolist()
    #return [list(x) for x in map(None, *M)]

def test():
    # Need to test both the C and python versions.
    import numpy
    import operator
    
    X = [[0, 1, 2], [4, 0, 5]]
    
    tests = [
        # 0
        (dim, (X,), {}, (2, 3)),                      # dim method
        (log, (10,), {"base":10}, 1.0),               # log_item
        (log, ([1,2,8],),{"base":2},[0.0,1.0,3.0]),   # log_list
        (log, (0,), {"safe":1}, 0.0),                 # log safe
        (mean, (1,), {}, "exception"),                # error on mean_item
        # 5
        (mean, (X[0],), {}, 1.0),                     # mean_list
        (mean, (X,), {}, [1.0, 3.0]),                 # mean_matrix, byrow
        (mean, (X,), {"byrow":0}, [2.0,0.5,3.5]),     # mean_matrix, bycol
        (median, (1,), {}, "exception"),              # error on median_item
        (median, (X[0],), {}, 1),                     # median_list
        # 10
        (median, ([2, 1, 3, 4],), {}, 2.5),           # median, unsorted list
        (median, (X,), {}, [1, 4]),                   # median_matrix, byrow
        (median, (X,), {"byrow":0}, [2.0,0.5,3.5]),   # median_matrix, bycol
        (var, (1,), {}, "exception"),                 # error on var_item
        (var, (X[0],), {}, 1.0),                      # var_list
        # 15
        (var, (X,), {}, [1.0, 7.0]),                  # var_matrix, byrow
        (var, (X,), {"byrow":0}, [8.0, 0.5, 4.5]),    # var_matrix, bycol
        (stddev, (1,), {}, "exception"),              # error on stddev_item
        (stddev, (X[0],), {}, 1.0),                   # stddev_list
        (stddev, (X,), {}, [1.0, 2.6]),               # stddev_matrix, byrow
        # 20
        (stddev, (X,), {"byrow":0}, [2.8, 0.7, 2.1]), # stddev_matrix, bycol
        (cov, (0,), {}, "exception"),                 # error on cov_item
        (cov, (X[0], X[1]), {}, 0.5),                 # cov_list
        (cov, (X,), {"byrow":1}, [[1.0,0.5],[0.5,7.0]]),  # cov_matrix, byrow
        # NotImplemented, cov_matrix, byrow=0
        (cor, (0,), {}, "exception"),                 # error on cor_item
        # 25
        (cor, (X[0], X[1]), {}, 0.2),                 # cor_list
        (cor, (X,), {"byrow":1}, [[1.0,0.2],[0.2,1.0]]),  # cor_matrix, byrow
        (cor, ([0,0], [1,2]), {}, "exception"),       # cor with 0 var
        (cor, ([0,0], [1,2]), {"safe":1}, 0.0),       # cor safe=0
        (cor, ([0,5,4], [3,2,1]), {"abs":0}, -0.8),   # cor abs=0
        # 30
        (cor, ([0,5,4], [3,2,1]), {"abs":1}, 0.8),    # cor abs=1
        (order, (0,), {}, "exception"),                 # error on order_item
        (order, (X[1],), {}, [1, 0, 2]),                # order_list
        (order, (X,), {"byrow":1}, [[0,1,2],[1,0,2]]),  # order_matrix, byrow
        (order, (X,0), {}, [[0,1,0],[1,0,1]]),          # order_matrix, bycol
        # 35
        (order, ([2, 3, 1],), {}, [2, 0, 1]),           # order vs rank
        (rank, ([2, 3, 1],), {}, [1.0, 2.0, 0.0]),
        (rank, ([1,2,1,1],), {}, [1.0,3.0,1.0,1.0]),    # rank with ties
        (rank, ([1,1,1,1],), {}, [1.5,1.5,1.5,1.5]),    # rank with ties
        (rank, ([1,2,1,1,2],), {}, [1.0,3.5,1.0,1.0,3.5]), # rank with ties
        # 40
        (rank, (X,), {}, [[0.0,1.0,2.0],[1.0,0.0,2.0]]),   # rank byrow=1
        (rank, (X,0), {}, [[0.0,1.0,0.0],[1.0,0.0,1.0]]),  # rank byrow=0
        (fisher_z, (0.2, 10), {}, 0.5),                 # fisher z
        (fisher_z, (0.2, 20), {}, 0.8),
        (fisher_z, (-0.5, 10), {}, -1.5),
        # 45
        (fisher_z, ([0.2, -0.5],10), {}, [0.5, -1.5]),
        (fisher_z, (0, 0), {}, "exception"),            # fisher z, N=0
        (transpose, (X,), {}, [[0,4],[1,0],[2,5]]),     # transpose
        (flatten, (X[0],), {}, [0, 1, 2]),              # flatten list
        (flatten, (X,), {}, [0, 1, 2, 4, 0, 5]),        # flatten matrix
        # 50
        (norm_mv, ([1, 2], 1, 2), {}, [0.0, 2.0]),
        (norm_mv, (X,1,2), {"byrow":1}, [[-0.4,1.0,2.4],[1.5,-0.6,2.1]]),
        (norm_mv, (X,1,2), {"byrow":0}, [[0.0,2.0,0.0],[2.0,0.0,2.0]]),
        (dnorm, (0, 0, 1), {}, 0.4),
        (dnorm, (0, 0, 2), {}, 0.2),
        # 55
        (dnorm, (5, 2, 6), {}, 0.1),
        (dnorm, (-1, 0, 2), {}, 0.2),
        (dnorm, ([0, 1, 2, 3], 0, 1), {}, [0.4,0.2,0.1,0.0]),   # dnorm_list
        (dnorm, (0, 0, 0), {}, "exception"),
        (max, (X,), {}, [2, 5]),                                # max
        # 60
        (max, (X,), {"byrow":0}, [4, 1, 5])                     # max, bycol
        ]

    def _fmt_single(x):
        if type(x) in [type(0.0), numpy.float64]:
            return float("%.1f" % x)
        return x
    def _fmt_list(X):
        x = [_fmt_single(x) for x in X]
        if type(X) is type(()):
            x = tuple(x)
        return x
    def _fmt_matrix(X):
        x = [_fmt_list(x) for x in X]
        if type(X) is type(()):
            x = tuple(x)
        return x
            
    for zzz, (fn, args, keywds, gold_standard) in enumerate(tests):
        #if zzz not in [33, 34, 35, 36]:
        #    continue
        try:
            output = fn(*args, **keywds)
        except Exception, x:
            output = "exception"
            if output != gold_standard:
                raise
            #raise
        status = "PASSED"

        output = _dispatch(output, _fmt_single, _fmt_list, _fmt_matrix)
        if str(output) != str(gold_standard):
            status = "FAILED"
        x = zzz, args, gold_standard, output, status
        print "\t".join(map(str, x))


# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    # This line doesn't work for some reason.  Will still use Python
    # versions.
    #from cjmath import *
    #raise ImportError
    import cjmath
except ImportError:
    pass
else:
    import sys
    this_module = sys.modules[__name__]
    this_dict, cjmath_dict = this_module.__dict__, cjmath.__dict__
    for name in cjmath_dict:
        if name.startswith("__"):
            continue
        if name in this_dict:
            this_dict["py_"+name] = this_dict[name]
        this_dict[name] = cjmath_dict[name]

if __name__ == '__main__':
    test()
