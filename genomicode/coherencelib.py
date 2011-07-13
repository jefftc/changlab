"""

Functions:
filter_genes
uniquify_genes
apply_standard_filters
apply_standard_filters_I

calc_coherence_score
calc_coherence_score_from_matrix
calc_coherence_p
calc_coherence_p_from_geo
load_alpha_distribution
estimate_alpha_distribution
choose_z_cutoff
fit_z_cutoff_to_ALPHA


Classes:
CoherenceScore
AlphaDistribution

"""
import os, sys

import arrayio
import geo_format
if geo_format not in arrayio.FORMATS:
    arrayio.FORMATS.append(geo_format)

if sys.hexversion < 0x02040000:
    def sorted(X, reverse=False):
        if type(X) is type({}):
            X = X.keys()
        X = X[:]
        X.sort()
        if reverse:
            X.reverse()
        return X

class CoherenceScore:
    """
    NROW
    NCOL
    CORS       Matrix of correlations, symmetric.
    ZCORS      Matrix of Fisher's Z scores, symmetric.
    PCORS      Number of positively correlated genes.
    NCORS      Number of negatively correlated genes.
    ALPHA      Number of positively correlated genes.
    MEDOID_I   Index of medoid, relative to CORS.
    PCORS_I    Index of PCOR genes, relative to CORS.
    NCORS_I    Index of NCOR genes, relative to CORS.
    
    """
    def __init__(self, NROW, NCOL, CORS, ZCORS, MEDOID_I, PCORS_I, NCORS_I):
        self.NROW, self.NCOL = NROW, NCOL
        self.CORS = CORS
        self.ZCORS = ZCORS
        self.PCORS = len(PCORS_I)
        self.NCORS = len(NCORS_I)
        #self.ALPHA = self.PCORS + self.NCORS
        self.ALPHA = self.PCORS
        self.MEDOID_I = MEDOID_I
        self.PCORS_I = PCORS_I
        self.NCORS_I = NCORS_I

class AlphaDistribution:
    def __init__(self, geneset_size, num_samples, alphas, counts):
        assert num_samples == sum(counts)
        self.geneset_size = geneset_size
        self.num_samples = num_samples
        # Sort alphas in increasing order.
        schwartz = zip(alphas, counts)
        schwartz.sort()
        self.alphas = [x[0] for x in schwartz]
        self.counts = [x[1] for x in schwartz]

def _find_medoid(CORS, CUTOFF):
    import jmath

    nrow, ncol = jmath.nrow(CORS), jmath.ncol(CORS)
    assert nrow == ncol
    
    pcors = [0]*nrow
    ncors = [0]*nrow
    for i in range(1, nrow):
        for j, x in enumerate(CORS[i][i+1:]):
            j += i
            if x >= CUTOFF:
                pcors[i] += 1
                pcors[j] += 1
            elif x <= -CUTOFF:
                ncors[i] += 1
                ncors[j] += 1

    max_cors = max_i = None
    for i in range(len(pcors)):
        num_cors = max(pcors[i], ncors[i])
        if max_cors is None or num_cors > max_cors:
            max_cors, max_i = num_cors, i
    return max_i

def _estimate_R_from_z(Z_SCORE, N, DELTA=0.001):
    # Given a Z-score, estimate the R value.
    import jmath
    
    i = 0
    min_R, max_R = -1.0, 1.0
    while min_R < max_R:
        R = (min_R + max_R) / 2.0
        z = jmath.fisher_z_item(R, N)
        if abs(z-Z_SCORE) < DELTA:
            break
        elif z > Z_SCORE:
            max_R = R - 1E-5
        elif z < Z_SCORE:
            min_R = R + 1E-5
        else:
            raise AssertionError, "How did I get here?"
        i += 1
    return R

def _cor_py(X, byrow=1, safe=1):
    import math
    
    assert X
    nrow, ncol = len(X), len(X[0])
    assert nrow and ncol
    if not byrow:
        raise NotImplementedError

    # Rounding error.  See below.
    EPS = 1E-5

    # If only 1 column, should set correlations to default of 0, or
    # there will be a divide-by-zero error.
    if ncol == 1:
        if not safe:
            raise AssertionError, "only 1 column in correlation."
        CORS = [[0]*nrow for i in range(nrow)]
        return CORS

    # Calculate the means and standard deviations of each row.
    for i in range(nrow):
        mean = float(sum(X[i])) / ncol

        # Covariance = E((X-mean(X))(Y-mean(Y)))
        # Subtract out the mean of each row for: X-mean(X)
        X[i] = [x-mean for x in X[i]]
        s = sum([x*x for x in X[i]])
        sd = math.sqrt(float(s)/(ncol-1))
        X[i] = [x/sd for x in X[i]]

    # Calculate the correlations.
    CORS = [[None]*nrow for i in range(nrow)]
    for i1 in range(nrow):
        X1 = X[i1]
        for i2 in range(i1, nrow):
            X2 = X[i2]
            s = 0
            for j in range(ncol):
                s += X1[j]*X2[j]
            x = float(s)/(ncol-1)
            
	    # Some values might be greater than 1 or -1 due to
	    # rounding error.  Fix this.
            if x >= 1.0 and x < 1.0+EPS:
                x = 1.0
            elif x < -1.0 and x > -1-EPS:
                x = -1.0
            CORS[i1][i2] = x
            CORS[i2][i1] = x
    return CORS

def _cor_py2(X, threshold, byrow=1, safe=1):
    import math
    
    assert X
    nrow, ncol = len(X), len(X[0])
    assert nrow and ncol
    if not byrow:
        raise NotImplementedError

    # If only 1 column, should set correlations to default of 0, or
    # there will be a divide-by-zero error.
    if ncol == 1:
        if not safe:
            raise AssertionError, "only 1 column in correlation."
        CORS = [[0]*nrow for i in range(nrow)]
        return CORS

    # Covariance = E((X-mean(X))(Y-mean(Y)))
    # Correlation = Cov(X, Y)/(SD(X) * SD(Y))

    # Calculate the means and standard deviations of each row.
    for i in range(nrow):
        mean = float(sum(X[i])) / ncol
        s = sum([(x-mean)*(x-mean) for x in X[i]])
        sd = math.sqrt(float(s)/(ncol-1))

        # Normalize each row to N(0, 1) to make the the correlation
        # coefficient calculation easier.
        X[i] = [(x-mean)/sd for x in X[i]]

    #X_large = [[None]*len(x) for x in X]
    #for i in range(len(X)):
    #    for j in range(len(X[i])):
    #        x = [abs(x) for x in X[i][j:]]
    #        X_large[i][j] = max(x)
    #X_large2 = [None] * ncol
    #for i in range(ncol):
    #    x = max([x[i] for x in X_large])
    #    X_large2[i] = x*x * (ncol-i)
        
    #X_sign = [None] * len(X)
    #for i in range(len(X)):
    #    X_sign[i] = [int(x >= 0) for x in X[i]]

    # Calculate the correlations.
    CORS = [[None]*nrow for i in range(nrow)]
    for i1 in range(nrow):
        X1 = X[i1]
        for i2 in range(i1, nrow):
            X2 = X[i2]

            #for j in range(ncol):
            #    if 0:
            #        pass
            t_12 = float(threshold)*(ncol-1)
            s = 0.0
            for j in range(ncol):
                #if s < t_12 - X_large[i1][j]*X_large[i2][j]*(ncol-j):
                #if s < t_12 - X_large2[j]:
                #    print "BREAK", j, ncol
                    #print "BREAK", i1, i2, s, t_12, X_l[i1][j], X_l[i2][j], \
                    #      j, ncol
                    #print X1[j], X2[j]
                    #break
                s += X1[j]*X2[j]
                #if s >= t_12:
                #    break
                
            #if i1 == 0 and i2 == 22:
            #    print "HERE", t_12, s, float(s)/(ncol-1)
            #    sys.exit(0)

            x = float(s >= t_12)
            CORS[i1][i2] = x
            CORS[i2][i1] = x
    return CORS

def cor(X, byrow=1, safe=1):
    import jmath
    import cjmath

    #return jmath.cor(X, byrow=byrow, safe=safe)   # numpy

    CORS1 = jmath.cor(X, byrow=byrow, safe=safe)   # numpy
    CORS2 = cjmath.cor_matrix(X, byrow=byrow, safe=safe)  # C or intrin
    CORS3 = _cor_py(X, byrow=byrow, safe=safe)     # pure python
    #print CORS1[0][:10]
    #print CORS2[0][:10]
    #assert jmath.equal_matrix(CORS1, CORS2)
    print jmath.assert_equal_matrix(CORS1, CORS2)
    print jmath.assert_equal_matrix(CORS1, CORS3)
    print jmath.assert_equal_matrix(CORS2, CORS3)
    print "EQUAL!"

    #threshold = 0.5
    #CORS3 = [None] * len(CORS1)
    #for i in range(len(CORS1)):
    #    CORS3[i] = [float(x >= threshold) for x in CORS1[i]]
    #CORS4 = _cor_py2(X, threshold, byrow=byrow, safe=safe)
    ##print CORS3[0][:10]
    ##print CORS4[0][:10]
    ##print CORS1[0][22], CORS2[0][22]
    #jmath.assert_equal_matrix(CORS3, CORS4)
    #print "EQUAL 34!"
    
    return CORS1


def _calc_coherence_score_h(X, Z_CUTOFF):
    # X is gene x sample matrix.
    import jmath

    #return None, None, [0], [0], [0]   # DEBUG
    
    # Need at least 2 samples, or no correlations.
    assert jmath.ncol(X) >= 2

    # Find the R-value that corresponds to Z_CUTOFF.
    R_CUTOFF = _estimate_R_from_z(Z_CUTOFF, jmath.ncol(X))

    # Find the medoid gene (correlated with the most neighbors).
    CORS = jmath.cor(X, byrow=1, safe=1)
    #CORS = cor(X, byrow=1, safe=1)
    MEDOID_I = _find_medoid(CORS, R_CUTOFF)

    # Find genes correlated with the medoid.
    cors_med = CORS[MEDOID_I]
    PCORS_I = [i for i in range(len(cors_med)) if cors_med[i] >= R_CUTOFF]
    NCORS_I = [i for i in range(len(cors_med)) if cors_med[i] <= -R_CUTOFF]
    # Arbitrarily set positive correlation as the group with the most
    # genes.
    if len(PCORS_I) < len(NCORS_I):
        PCORS_I, NCORS_I = NCORS_I, PCORS_I

    if 0:
        # DEBUGGING
        import ccoherence
        x = ccoherence._calc_coherence_score(X, Z_CUTOFF)
        (C_CORS, C_ZCORS, C_MEDOID_I,
         C_ALPHA_I, C_BETA_I, C_len_ALPHA_I, C_len_BETA_I) = x
        print "C"
        print len(C_CORS), len(C_CORS[0])
        print C_CORS[0][:10]
        print C_CORS[MEDOID_I][:10]
        print C_MEDOID_I, C_ALPHA_I, C_BETA_I
        print ccoherence._find_medoid(C_CORS, R_CUTOFF)
        print "PYTHON"
        print len(CORS), len(CORS[0])
        print CORS[0][:10]
        print CORS[MEDOID_I][:10]
        print MEDOID_I, ALPHA_I, BETA_I
        print ccoherence._find_medoid(CORS, R_CUTOFF)
        print "PY R_CUTOFF IS", R_CUTOFF
        sys.exit(0)

    return X.nrow(), X.ncol(), CORS, None, MEDOID_I, PCORS_I, NCORS_I

def calc_coherence_score(dataset, geneset, z_cutoff):
    # Return a CoherenceScore object.
    import jmath
    I = geneset.get_indexes(dataset)
    # This can happen if the gene_id_name of the GeneSet doesn't match
    # the gene IDs of the data set.
    assert I, "Dataset contains no genes from geneset."
    X = dataset.slice(row=I)
    assert jmath.nrow(X) > 0 and jmath.ncol(X) > 0, "Empty matrix"
    assert jmath.ncol(X) >= 2, "not enough samples (%d)" % jmath.ncol(X)
    return calc_coherence_score_from_matrix(X, z_cutoff)

def calc_coherence_score_from_matrix(X, z_cutoff):
    # Return a CoherenceScore object.
    x = _calc_coherence_score_h(X, z_cutoff)
    return CoherenceScore(*x)

def _linreg(X, Y):
    """
    Summary
    Linear regression of y = ax + b
    Usage
    real, real, real = linreg(list, list)
    Returns coefficients to the regression line "y=ax+b" from x[] and
    y[], and R^2 Value
    """
    from math import sqrt
    assert len(X) == len(Y)
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in zip(X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    for x, y in zip(X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    RR = 1.0 - residual/meanerror
    #ss = residual / (N-2)
    #Var_a, Var_b = ss * N / det, ss * Sxx / det
    return a, b, RR

def _calc_p_from_samples(obs, samples):
    num_samples = len(samples)
    #multiplier = 10**max((len(str(num_samples))-1), 2)
    
    # P-value is probability of sampling an observation at least as
    # extreme as obs.
    num_extreme = len([x for x in samples if x >= obs])
    mod = "="
    if num_extreme == 0:
        # No 0 p-value, and valid number for nl10p.
        num_extreme = 1
        mod = "<"
    pvalue = float(num_extreme) / num_samples
    return num_extreme, num_samples, mod, pvalue

def calc_coherence_p(dataset, geneset, z_cutoff, obs_alpha, num_samples=100):
    alphas = [None] * num_samples
    for i in range(num_samples):
        x = calc_coherence_score(dataset, geneset, z_cutoff)
        alphas[i] = x.ALPHA
    x = _calc_p_from_samples(obs_alpha, alphas)
    num_extreme, num_samples, mod, pvalue = x
    return mod, pvalue

def calc_coherence_p_from_dist(obs_alpha, num_probes, size2dist, sizes=None):
    # Can set the sizes that should be used for interpolation (for
    # debugging).
    # Interpolate between two sizes that have been sampled.  Will not
    # extrapolate outwards.
    import jmath
    parzen_h = 0.5

    assert num_probes > 0, "must have at least 1 gene"
    assert obs_alpha >= 1 and obs_alpha <= num_probes
    
    # If the user specific specific sizes to use, then make sure they
    # have been sampled.
    if sizes is not None:
        for s in sizes:
            assert s in size2dist

    # If this size is already sampled, then no need to interpolate.
    #if num_probes in size2dist and (sizes is None or num_probes in sizes):
    if sizes is None and num_probes in size2dist:
        sizes = [num_probes]

    # Select the sizes to use to interpolate num_probes.
    if sizes is None:
        # Strategy: Interpolate using the next larger and next smaller
        # size.
        smaller = [s for s in size2dist if s <= num_probes]
        larger = [s for s in size2dist if s > num_probes]
        assert smaller, \
               "Gene set [%d] is too small for p-value estimation [%d]." % (
            num_probes, min(size2dist))
        assert larger, \
               "Gene set [%d] is too big for p-value estimation [%d]." % (
            num_probes, max(size2dist))
        smaller, larger = max(smaller), min(larger)
        sizes = [smaller, larger]
        ## Strategy: Use the two signature sizes that are closest.
        ## BUG: Should use the closest smaller and closest larger one.
        #schwartz = [(abs(s-num_probes), s) for s in size2dist]
        #schwartz.sort()
        #sizes = [x[-1] for x in schwartz]
        #sizes = sizes[:2]

    if len(sizes) == 1:
        # Only one size, must be already sampled.
        size = sizes[0]
        assert size in size2dist
        dist = size2dist[size]
        N = dist.num_samples
        eps = 1.0 / N / 2
        pvalue = 1.0 - jmath.pparzen(
            obs_alpha, dist.alphas, dist.counts, h=parzen_h, eps=eps)
    elif 0:
        # Estimate the P-value of ALPHA from the data in sizes.  Try
        # different p-values in a binary search strategy.  For each
        # p-value:
        # 1.  Find the ALPHA at each gene set size.
        # 2.  Find ALPHA for this gene set size (using regression).
        # 3.  Compare against observed ALPHA.
        N = min([size2dist[s].num_samples for s in sizes])
        PVALUE_DELTA = 1.0/N/2
        # Search for good bounds for the p-value by looking at the
        # other sizes.  Doesn't work: range is [0, 1].
        #pvalues = []
        #for dist in size2dist.values():
        #    N = dist.num_samples
        #    eps = 1.0 / N / 2
        #    pvalue = 1.0 - jmath.pparzen(
        #        obs_alpha, dist.alphas, dist.counts, h=parzen_h, eps=eps)
        #    pvalues.append(pvalue)
        #print min(pvalues), max(pvalues)
        pvalue_min = 0
        pvalue_max = 1
        while pvalue_min < pvalue_max:
        #while abs(pvalue_min-pvalue_max) > PVALUE_DELTA:  # slower
            pvalue = (pvalue_min + pvalue_max) / 2.0

            x = _fit_size_to_alpha(pvalue, size2dist, sizes, parzen_h=parzen_h)
            x0, x1 = x

            # Interpolate the alpha at geneset.size.
            p_alpha = x0 + x1 * num_probes
            #print pvalue, pvalue_min, pvalue_max, p_alpha, obs_alpha

            if p_alpha < obs_alpha:
                # Predicted alpha is lower than desired one, so lower the
                # p-value.
                pvalue_max = pvalue - PVALUE_DELTA
            else:
                pvalue_min = pvalue + PVALUE_DELTA
    else:
        N = min([size2dist[s].num_samples for s in sizes])
        pvalue = _find_p_value(
            obs_alpha, num_probes, size2dist, sizes=sizes, parzen_h=parzen_h)

    mod = "="
    if pvalue < 1.0/N:
        # Not sampled enough to estimate p-value.
        mod = "<"
    # Round the p-value to the proper resolution.
    pvalue = max(round(pvalue*N), 1.0)/N
    return mod, pvalue

def calc_coherence_p_from_geo(
    geo_path, GSEID, GPLID, num_probes, obs_alpha, sizes=None):
    # Load the alpha distribution file.
    size2dist = load_alpha_distribution(geo_path, GSEID, GPLID)
    x = calc_coherence_p_from_dist(
        obs_alpha, num_probes, size2dist, sizes=sizes)
    return x

def _find_alpha_range(X_obs, X_count, h=1, eps=1E-10):
    import math
    
    N = sum(X_count)
    A = eps*N*h*math.sqrt(2*math.pi)/max(X_count)
    B = h*math.sqrt(-2*math.log(A))
    x_min, x_max = min(X_obs)-B, max(X_obs)+B
    return x_min, x_max

def _increment_alpha(x_cur, x_max, q_cur, q, X_obs, X_count, h=1, eps=1E-10):
    import jmath
    #x_delta = (x_max-x_cur)/2.0

    #print "INC ALPHA"
    #num_failures = 0
    #while abs(q_cur-q) > eps:
    x = x_cur
    while q-q_cur > eps:
        #x = x_cur + x_delta
        #if x >= x_max:
        #    x = x_max
        #    break
        x = (x_cur + x_max)/2.0
        d = jmath._pparzen_ul(x_cur, x, X_obs, X_count, h=h, eps=eps)
        if d < eps:
            break
        #print q_cur, q, d, x_cur, x, x_max, x_delta
        #print q_cur, q, d, x_cur, x, x_max
        if q_cur + d < q:
            q_cur += d
            x_cur = x
            #num_failures = 0
        else:
            #x_delta = x_delta / (2.0+num_failures/2.0)
            x_max = x
            #num_failures += 1
    return x

def _find_p_value(
    alpha_obs, size_obs, size2dist, sizes=None, parzen_h=0.5):
    if sizes is None:
        sizes = size2dist.keys()

    # Figure out the right eps value.
    N = min([x.num_samples for x in size2dist.values()])
    eps = 1.0 / N / 2.0

    # For each size, calculate the alpha that corresponds to q=0.
    len_sizes = len(sizes)
    alphas = [None] * len_sizes
    max_alphas = [None] * len_sizes
    for i, size in enumerate(sizes):
        dist = size2dist[size]
        x = _find_alpha_range(dist.alphas, dist.counts, h=parzen_h, eps=eps)
        alpha_min, alpha_max = x
        alphas[i] = alpha_min
        max_alphas[i] = alpha_max

    q_min, q_max = 0.0, 1.0
    while q_max-q_min > eps:
        q = (q_min+q_max)/2.0
        
        # For each size, increment the alpha until the quantile reaches q.
        alphas_new = [None] * len(sizes)
        for i in range(len_sizes):
            dist = size2dist[sizes[i]]
            x = _increment_alpha(
                alphas[i], max_alphas[i], q_min, q, dist.alphas, dist.counts,
                h=parzen_h, eps=eps)
            alphas_new[i] = x

        x1, x0, RR = _linreg(sizes, alphas_new)
        alpha = x0 + x1 * size_obs
        if alpha < alpha_obs:
            q_min = q
            alphas = alphas_new
        elif alpha > alpha_obs:
            q_max = q
        else:
            break
    pvalue = 1.0 - q
    return pvalue
    
def _fit_size_to_alpha(pvalue, size2dist, sizes=None, parzen_h=0.5):
    # Do a regression to fit: A = x0 + x1*S
    # Return x0, x1.
    import jmath
    
    if sizes is None:
        sizes = size2dist.keys()
    q = 1.0 - pvalue
    alphas = []
    for size in sizes:
        dist = size2dist[size]
        eps = 1.0 / dist.num_samples / 2
        a = jmath.qparzen(q, dist.alphas, dist.counts, h=parzen_h, eps=eps)
        alphas.append(a)
    x1, x0, RR = _linreg(sizes, alphas)
    return x0, x1

def load_alpha_distribution(geo_path, GSEID, GPLID):
    # Return dictionary of geneset size -> AlphaDistribution
    import geolib
    import iolib

    filename = geolib.FileFactory.find_alpha_dist_file(geo_path, GSEID, GPLID)
    assert filename, "Missing alpha dist file for %s-%s" % (GSEID, GPLID)
    x = iolib.split_tdf(open(filename).read())
    
    #x = [map(int, x) for x in x]
    # Using GenericObjects is very expensive.
    #data = [filelib.GenericObject(
    #    size=x[0], samples=x[3], count=x[2], alpha=x[1]) for x in x]
    #for d in data:
    #    size2data.setdefault(d.size, []).append(d)

    if sys.hexversion >= 0x02040000:
        # This version is faster, but only available on 2.4 and above.
        import itertools
        size2data = itertools.groupby(x, lambda x: x[0])
        iter_size2data = size2data
    else:
        size2data = {}
        for x in x:
            # size, alpha, count, samples
            #x = map(int, x)
            size = x[0]
            if size not in size2data:
                size2data[size] = []
            size2data[size].append(x)
        iter_size2data = size2data.iteritems()
        
    size2dist = {}
    for size, ds in iter_size2data:
        ds = list(ds)
        #samples = ds[0].samples
        samples = ds[0][3]
        for d in ds:
            assert d[3] == samples
            #assert d.samples == samples
        #alphas = [d.alpha for d in ds]
        #counts = [d.count for d in ds]
        alphas = [int(x[1]) for x in ds]
        counts = [int(x[2]) for x in ds]
        size = int(size)
        x = AlphaDistribution(size, int(samples), alphas, counts)
        size2dist[size] = x
    return size2dist

def estimate_alpha_distribution(size2dist, num_probes, sizes=None):
    parzen_h = 0.5

    assert num_probes > 0, "must have at least 1 gene"

    # Load the alpha distribution file.
    #size2dist = load_alpha_distribution(geo_path, GSEID, GPLID)

    # If the user specific specific sizes to use, then make sure they
    # have been sampled.
    if sizes is not None:
        for s in sizes:
            assert s in size2dist

    # If this size is already sampled, then no need to interpolate.
    if sizes is None and num_probes in size2dist:
        return size2dist[num_probes]

    # Select the sizes to use to interpolate num_probes.
    if sizes is None:
        # Strategy: Interpolate using the next larger and next smaller
        # size.
        smaller = [s for s in size2dist if s <= num_probes]
        larger = [s for s in size2dist if s > num_probes]
        assert smaller, \
               "Gene set [%d] is too small for estimation [%d]." % (
            num_probes, min(size2dist))
        assert larger, \
               "Gene set [%d] is too big for estimation [%d]." % (
            num_probes, max(size2dist))
        smaller, larger = max(smaller), min(larger)
        sizes = [smaller, larger]

    assert len(sizes) >= 2, "Must interpolate from at least two sizes."
    
    # Estimate the count of ALPHA from the data in sizes.
    big_dist = size2dist[max(sizes)]
    small_dist = size2dist[min(sizes)]
    alphas = range(min(small_dist.alphas), max(big_dist.alphas))
    # pvalue is tuple of ('=', 0.7027).  Just want the number.
    x = [calc_coherence_p_from_dist(x, num_probes, size2dist)
         for x in alphas]
    pvalues = [x[1] for x in x]
    # Convert the p-values to counts.
    TOTAL_COUNTS = 1000000
    counts = [None] * len(pvalues)
    counts[0] = int(round((1.0-pvalues[0])*TOTAL_COUNTS))
    for i in range(1, len(pvalues)):
        counts[i] = int(round((pvalues[i-1]-pvalues[i])*TOTAL_COUNTS))
    dist = AlphaDistribution(num_probes, sum(counts), alphas, counts)
    return dist

def choose_z_cutoff(dataset, gene_id_name=None, debug_handle=None):
    # Choose a reasonable Z cutoff for this data set.  Z-cutoff is the
    # 95% Z-score of a set of random genes.  (A correlation
    # coefficient that exceeds 95% of the correlations in this gene
    # set.)
    import arrayio
    import jmath
    import GeneSet

    QUANTILE = 0.95
    NUM_GENES = 1000
    
    gene_id_name = gene_id_name or arrayio.GENE_ID

    geneset = GeneSet.RandomGeneSet(gene_id_name, dataset, NUM_GENES)
    X = dataset.slice(row=geneset.get_genes(),  row_header=gene_id_name)
    N = jmath.ncol(X)
    assert N > 2, "insufficient samples"

    ACORS = jmath.cor(X, byrow=1, safe=1, abs=1)
    cors_s = []
    for i in range(jmath.nrow(ACORS)-1):
        x = ACORS[i][i+1:]
        cors_s.extend(x)
    TOTAL = (jmath.nrow(ACORS)*(jmath.nrow(ACORS)-1))/2
    assert len(cors_s) == TOTAL, "%d %d %d" % (
        jmath.nrow(ACORS), len(cors_s), TOTAL)
    cors_s.sort()

    i = int(round(len(cors_s) * QUANTILE))
    R = cors_s[i]

    Z_CUTOFF = jmath.fisher_z_item(R, N)
    return Z_CUTOFF

def fit_z_cutoff_to_ALPHA(dataset, geneset, ALPHA, debug_handle=None):
    # Find a Z cutoff such that a data set will get a specific ALPHA
    # score for a gene set.
    assert ALPHA >= 1 and ALPHA <= len(geneset)
    DELTA = 1E-5
    Z_min, Z_max = DELTA, 50.0
    i = 0
    found = 0
    while Z_min < Z_max:
        Z = (Z_min + Z_max) / 2.0
        cscore = calc_coherence_score(dataset, geneset, Z)
        if debug_handle:
            x = i, Z_min, Z_max, Z, cscore.ALPHA, ALPHA
            print >>debug_handle, "\t".join(map(str, x))
        if abs(cscore.ALPHA-ALPHA) < DELTA:
            found = 1
            break
        elif cscore.ALPHA > ALPHA:
            # If this ALPHA is too large, then I need to increase the
            # stringency by increasing the Z cutoff.
            Z_min = Z + DELTA
        elif cscore.ALPHA < ALPHA:
            Z_max = Z - DELTA
        else:
            raise AssertionError, "How did I get here?"
        i += 1
    assert found, "I could not fit a good Z cutoff."
    return Z


def test_calc_coherence_score():
    import matrixfns as mf
    import ccoherence
    import numpy
    import time

    X = [1, 2, 3, 4, 5, 9, 15]
    print ccoherence.mean(X)   #  5.57
    print ccoherence.var(X)    # 23.95

    
    Z_CUTOFF = 4
    X = [[1, 2, 3], [2, 5, 4], [1, 7, 8], [9, 5, 7]]
    #X = "hello"

    print "TESTING COV"
    print numpy.cov(X, None, 1, 0)
    print ccoherence.cov(X, bycol=0)
    print
    
    print "TESTING COR"
    print numpy.corrcoef(X)
    print ccoherence.safe_cor(X, bycol=0)
    print
    
    cor = ccoherence.safe_cor(X, bycol=0)

    print "TESTING SCORE"
    print fisher_z_transform(cor, mf.ncol(X))
    print ccoherence.fisher_z_transform(cor, mf.ncol(X))
    print

    #print "TESTING SPLIT"
    #data = openfh("geo_db/annot/GPL96.txt.gz").read()
    #print [x.rstrip("\r\n").split("\t") for x in data.split("\n")][:5]
    #print ccoherence.split_tab_delimited_format(data)[:5]
    #print

    print "TIMING"
    start = time.time()
    for i in range(10000):
        numpy.corrcoef(X)
    print time.time()-start

    start = time.time()
    for i in range(10000):
        ccoherence.safe_cor(X, bycol=0)
    print time.time()-start

    start = time.time()
    for i in range(10000):
        fisher_z_transform(cor, mf.ncol(X))
    print time.time()-start
    
    start = time.time()
    for i in range(10000):
        ccoherence.fisher_z_transform(cor, mf.ncol(X))
    print time.time()-start

    #start = time.time()
    #for i in range(10):
    #    [x.rstrip("\r\n").split("\t") for x in data.split("\n")]
    #print time.time()-start
    
    #start = time.time()
    #for i in range(10):
    #    ccoherence.split_tab_delimited_format(data)
    #print time.time()-start

    # Check for memory leaks
    print "Running forever to check for memory leaks"
    for i in range(10000):
        for j in range(10000):
            calc_coherence_score()

def test_calc_coherence_p():
    geo_path = "geo_db"
    GSEID = "GSE3494"
    GPLID = "GPL96"
    len_geneset = 300
    alpha = 200

    x = _calc_coherence_p(geo_path, GSEID, GPLID, len_geneset, alpha)
    p_value, multiplier, mod, quant2alpha = x
    
    quants = quant2alpha.keys()
    quants.sort()
    for quant in quants:
        x = len_geneset, multiplier, quant, quant2alpha[quant]
        print "\t".join(map(str, x))

    x = GSEID, GPLID, len_geneset, alpha, multiplier, mod, p_value
    print "\t".join(map(str, x))

TDF_FILE = None
TDF_CONTENTS = None
def _read_tdf(filename):
    global TDF_FILE, TDF_CONTENTS
    import itertools
    import iolib

    if filename != TDF_FILE:
        x = [x for x in iolib.split_tdf(open(filename).read())]
        TDF_CONTENTS = list(itertools.chain(*x))
        TDF_FILE = filename
    return TDF_CONTENTS

def _get_unique_indexes(names, good_names):
    # Strategies to keep track of the indexes, in order of fastest to
    # slowest.
    # - Make a list of the maximum possible size, and keep track of
    #   the length of the list manually.
    # - Use the indexes saved in seen.values() and sort them.  (A tiny
    #   bit slower for 12,000 genes).
    # - Make a list of indexes and .append an index each time.  (A lot
    #   slower).
    indexes = [None] * len(names)
    num_indexes = 0
    seen = {}
    for i, name in enumerate(names):
        if name not in good_names:
            continue
        if name in seen:
            continue
        seen[name] = 1
        indexes[num_indexes] = i
        num_indexes += 1
    #indexes = sorted(seen.values())
    indexes = indexes[:num_indexes]
    return indexes

def apply_standard_filters_I(
    dataset, id_name, filename=None, allowed_gene_ids=None):
    # filename is a file containing a list of IDs to allow.
    # allowed_gene_ids is a list of the IDs to allow.  At least one of
    # these should be provided.
    import GeneSet
    
    assert filename or allowed_gene_ids
    
    gene_ids = allowed_gene_ids or []
    if filename is not None:
        gene_ids = _read_tdf(filename)

    # Make sure each gene_id is valid.
    #x1 = GeneSet.select_valid_gene_id_new(gene_ids)
    #x2 = GeneSet.select_valid_gene_id(gene_ids)
    #print "OLD", x1[:10]
    #print "NEW", x2[:10]
    #print x1 == x2
    gene_ids = GeneSet.select_valid_gene_id(gene_ids)
    assert gene_ids, "no genes remain after filter"
    
    # Turn gene_ids into a dictionary.
    # This line takes about 10% of the running time of this function.
    gene_ids = {}.fromkeys(gene_ids)

    # Figure out the indexes of the matrix that correspond to the
    # gene_ids.  Make sure there aren't any duplicates.
    names = dataset.row_names(id_name)
    indexes = _get_unique_indexes(names, gene_ids)
    
    ## seen = {}
    ## indexes = [None] * dataset.nrow()
    ## num_indexes = 0
    ## for i, id in enumerate(dataset.row_names(id_name)):
    ##     if id not in gene_ids:
    ##         continue
    ##     if id in seen:
    ##         continue
    ##     seen[id] = 1
    ##     indexes[num_indexes] = i
    ##     num_indexes += 1
    
    ## # Strategies to keep track of the indexes, in order of fastest to
    ## # slowest.
    ## # - Make a list of the maximum possible size, and keep track of
    ## #   the length of the list manually.
    ## # - Use the indexes saved in seen.values() and sort them.  (A tiny
    ## #   bit slower for 12,000 genes).
    ## # - Make a list of indexes and .append an index each time.  (A lot
    ## #   slower).
    ## #indexes = sorted(seen.values())
    ## indexes = indexes[:num_indexes]
    return indexes

def apply_standard_filters(
    dataset, id_name, filename=None, allowed_gene_ids=None):
    indexes = apply_standard_filters_I(
        dataset, id_name, filename=filename, allowed_gene_ids=allowed_gene_ids)
    return dataset.matrix(row=indexes)

## def apply_standard_filters(
##     dataset, id_name, filename=None, allowed_gene_ids=None):
##     # filename is a file containing a list of IDs to allow.
##     # allowed_gene_ids is a list of the IDs to allow.  At least one of
##     # these should be provided.
##     dataset = uniquify_genes(dataset, id_name)
##     dataset = filter_genes(
##         dataset, id_name, filename=filename, allowed_gene_ids=allowed_gene_ids)
##     return dataset

def filter_genes(dataset, id_name, filename=None, allowed_gene_ids=None):
    # filename is a file containing a list of IDs to allow.
    # allowed_gene_ids is a list of the IDs to allow.  At least one of
    # these should be provided.
    import itertools
    import GeneSet
    import iolib
    
    assert filename or allowed_gene_ids
    
    gene_ids = allowed_gene_ids or []
    if filename is not None:
        #for x in iolib.split_tdf(open(filename).read()):
        #    gene_ids.extend(x)
        x = [x for x in iolib.split_tdf(open(filename).read())]
        gene_ids = itertools.chain(*x)

    # Get a list of unique valid gene IDs.
    gene_ids = {}.fromkeys(gene_ids).keys()
    gene_ids = GeneSet.select_valid_gene_id(gene_ids)
    #gene_ids = [x for x in gene_ids if GeneSet.is_valid_gene_id(x)]
    assert gene_ids, "no genes remain after filter"

    return dataset.matrix(row=gene_ids, row_header=id_name)

def uniquify_genes(dataset, id_name):
    gene_ids = dataset.row_names(id_name)

    seen = {}
    indexes = [None] * len(gene_ids)
    num_indexes = 0
    for i, id in enumerate(gene_ids):
        if id in seen:
            continue
        seen[id] = i
        indexes[num_indexes] = i
        num_indexes += 1
    # Strategies to keep track of the indexes, in order of fastest to
    # slowest.
    # - Make a list of the maximum possible size, and keep track of
    #   the length of the list manually.
    # - Use the indexes saved in seen.values() and sort them.  (A tiny
    #   bit slower for 12,000 genes).
    # - Make a list of indexes and .append an index each time.  (A lot
    #   slower).
    #indexes = sorted(seen.values())
    indexes = indexes[:num_indexes]
    return dataset.matrix(row=indexes)




# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    # This line doesn't work for some reason.  Will still use Python
    # versions.
    #from ccoherence import *
    #raise ImportError
    import ccoherencelib
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in ccoherencelib.__dict__.keys():
        if not name.startswith("__"):
            this_module.__dict__[name] = ccoherencelib.__dict__[name]
