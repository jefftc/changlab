"""

Functions:
correct_background
normalize_within_arrays


These functions produce values that are slightly different than when
run in R.  It may be because of rounding errors when converting from R
variables to Python variables.

"""

def correct_background(
    red_signal, green_signal, red_back, green_back, method=None, offset=None):
    # Apply a background correction and return a tuple of (red_matrix,
    # green_matrix).  See the documentation for the backgroundCorrect
    # method in the R limma package for the allowable values for method
    # and offset.
    # 
    # Common ways to use this are:
    # method="normexp"   offset=50
    # method="subtract"
    import sys
    import StringIO
    import jmath

    # If method is incorrect, an exception will be raise in R.

    # Check the input matrices.
    assert red_signal
    n, m = len(red_signal), len(red_signal[0])
    matrices = [red_signal, green_signal, red_back, green_back]
    for mat in matrices:
        assert mat
        assert len(mat) == n
        for row in mat:
            assert len(row) == m

    R = jmath.start_R()
    # If limma library doesn't exist, an exception will be raised in
    # R.
    R('library("limma")')
    R('RG <- list()')
    jmath.R_equals(red_signal, "RG$R")
    jmath.R_equals(green_signal, "RG$G")
    jmath.R_equals(red_back, "RG$Rb")
    jmath.R_equals(green_back, "RG$Gb")
    #print jmath.R2py_matrix(R('RG$R[1:5,1:5]'))
    #print jmath.R2py_matrix(R('RG$Rb[1:5,1:5]'))
    args = ["RG"]
    if method:
        args.append('method="%s"' % method)
    if offset is not None:
        args.append("offset=%d" % offset)
    x = ", ".join(args)
    old_stdout = None
    try:
        old_stdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        # This will generate text to the console.  Capture it so it
        # doesn't mess up the screen.
        R('RG <- backgroundCorrect(%s)' % x)
    finally:
        if old_stdout is not None:
            sys.stdout = old_stdout
    R_matrix = jmath.R2py_matrix(R('RG$R'))
    G_matrix = jmath.R2py_matrix(R('RG$G'))
    return R_matrix, G_matrix

def normalize_within_arrays(red_signal, green_signal, method=None):
    # Normalize the signal within arrays and return a matrix of the
    # normalized signal values.  See the documentation for the
    # normalizeWithinArrays method in the R limma package for the
    # allowable values for method.
    # 
    # Common ways to use this are:
    # method="none"
    # method="loess"
    import jmath

    # If method is incorrect, an exception will be raise in R.

    # Check the input matrices.
    assert red_signal
    n, m = len(red_signal), len(red_signal[0])
    matrices = [red_signal, green_signal]
    for mat in matrices:
        assert mat
        assert len(mat) == n
        for row in mat:
            assert len(row) == m

    R = jmath.start_R()
    # If limma library doesn't exist, an exception will be raised in
    # R.
    R('library("limma")')
    R('RG <- list()')
    jmath.R_equals(red_signal, "RG$R")
    jmath.R_equals(green_signal, "RG$G")
    args = ["RG"]
    if method:
        args.append('method="%s"' % method)
    x = ", ".join(args)
    #print x
    R('MA <- normalizeWithinArrays(%s)' % x)
    x = jmath.R2py_matrix(R('as.matrix(MA)'))
    return x
