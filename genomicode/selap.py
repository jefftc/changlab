"""

Functions:
selap_make_raw
selap_predict_raw

find_selap

"""
import os, sys

def selap_make_raw(
    matrix_file, penalty, matlab_bin=None, selap_path=None, outpath=None):
    # matrix_file should contain a matrix (no headers) of samples x
    # pathway predictions.
    import matlab
    
    # Set defaults.
    matrix_file = os.path.realpath(matrix_file)
    assert os.path.exists(matrix_file)
    selap_path = find_selap(selap_path)
    assert selap_path is not None, "I could not find SELAP."
    selap_path = os.path.realpath(selap_path)
    outpath = outpath or "."

    lines = []
    w = lines.append
    w("addpath '%s/Utilities';\n" % selap_path)
    w("X = load('%s');\n" % matrix_file)
    w("X = logit(X);\n")
    w("[mu sig prob] = selapMixMap(X, %s);\n" % penalty)
    #w("predictions = mvnMixtureProb(X', mu, sig, prob);\n")
    w("save('mu.txt', 'mu', '-ASCII', '-TABS');\n")
    w("save('sig.txt', 'sig', '-ASCII', '-TABS');\n")
    w("save('prob.txt', 'prob', '-ASCII', '-TABS');\n")
    #w("save('predict.txt', 'predictions', '-ASCII', '-TABS');\n")
    script = "".join(lines)
    x = matlab.run(script, matlab_bin=matlab_bin, working_path=outpath)
    return x

def selap_predict_raw(
    matrix_file, mu_file, sig_file, prob_file, matlab_bin=None,
    selap_path=None, outpath=None):
    # matrix_file should contain a matrix (no headers) of samples x
    # pathway predictions.
    import matlab
    
    # Set defaults.
    matrix_file = os.path.realpath(matrix_file)
    mu_file = os.path.realpath(mu_file)
    sig_file = os.path.realpath(sig_file)
    prob_file = os.path.realpath(prob_file)
    
    assert os.path.exists(matrix_file)
    assert os.path.exists(mu_file)
    assert os.path.exists(sig_file)
    assert os.path.exists(prob_file)

    selap_path = find_selap(selap_path)
    assert selap_path is not None, "I could not find SELAP."
    selap_path = os.path.realpath(selap_path)
    outpath = outpath or "."

    lines = []
    w = lines.append
    w("addpath '%s/Utilities';\n" % selap_path)
    w("X = load('%s');\n" % matrix_file)
    w("X = logit(X);\n")

    w("mu = load('%s');\n" % mu_file)
    w("sig = load('%s');\n" % sig_file)
    w("prob = load('%s');\n" % prob_file)
    
    w("predictions = mvnMixtureProb(X', mu, sig, prob);\n")
    w("save('predict.txt', 'predictions', '-ASCII', '-TABS');\n")
    script = "".join(lines)
    x = matlab.run(script, matlab_bin=matlab_bin, working_path=outpath)
    return x

def find_selap(default_path):
    search_paths = [
        default_path,
        "/home/jchang/projects/genepattern/data/SELAPver3",
        "/Volumes/Users/jchang/projects/genepattern/data/SELAPver3",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        # Test for some BFRM files.
        files = [
            "runme.txt", "subanalysis.m", "Utilities",
            ]
        complete = True   # Is this distribution complete.
        for file in files:
            filename = os.path.join(spath, file)
            if not os.path.exists(filename):
                complete = False
                # Don't break here, so we can diagnose missing files.
                #break
        if not complete:
            continue
        path = spath
        break
    return path
