"""

Functions:
normalize
normalize_binreg    Requires binreg

"""
import os

def normalize(X, which_columns=None):
    # X is a Matrix of the data.  which_columns is a list of the
    # columns used to calculate the quantiles.  If None, then will use
    # every column.  Returns the normalized version of the matrix.
    import jmath

    # Make a list of the quantiles from smallest to largest.
    if which_columns is None:
        which_columns = range(X.ncol())
    X_quant = X[(None, which_columns)]
    quantiles = sorted(jmath.median(X_quant))

    X_norm = [[None]*X.ncol() for i in range(X.nrow())]
    for i in range(X.ncol()):
        # Normalize each column.
        col = X[(None, i)]
        O = jmath.order(col)
        for o, v in zip(O, quantiles):
            X_norm[o][i] = v

    Y = X.matrix()
    Y._X = X_norm
    return Y


def normalize_binreg(X, which_columns=None, matlab=None, binreg_path=None):
    # X is a Matrix of the data.  which_columns is a list of the
    # columns used to calculate the quantiles.  If None, then will use
    # every column.  Returns the normalized version of the matrix.
    #import copy
    import subprocess
    import tempfile
    
    import jmath
    import Matrix
    import binregfns
    from dwdnorm import _write_matlab_matrix, _parse_normalized_matrix
    from dwdnorm import _safe_unlink
    
    # Make a list of the quantiles from smallest to largest.
    if which_columns is None:
        which_columns = range(X.ncol())
    X_quant = X[(None, which_columns)]
    m = jmath.median(X_quant)
    m = Matrix.InMemoryMatrix([m])
    
    # Set defaults.
    matlab = matlab or "matlab"
    binreg_path = binregfns.find_binreg_20(binreg_path)
    assert binreg_path, "I could not find Binreg2.0"
    binreg_path = os.path.realpath(binreg_path)
    temp_path = "."

    # Start an instance of matlab.
    matlab_args = [
        "-nosplash", "-nodesktop", "-nodisplay", "-nojvm"]
    x = " ".join(matlab_args)
    cmd = "%s %s" % (matlab, x)
    #print cmd
    #w, r = os.popen4(cmd, bufsize=0)
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout

    # Run the normalization.
    X_file = m_file = None
    try:
        x, X_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        x, m_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        _write_matlab_matrix(X, X_file)
        _write_matlab_matrix(m, m_file)
        script = _format_exec_file(X_file, m_file, binreg_path)
        w.write(script)
        w.close()
        lines = r.readlines()
    finally:
        _safe_unlink(X_file)
        _safe_unlink(m_file)

    # Debug: print out the matrix returned by Matlab.
    #for l in lines:
    #    print l,
    
    X_norm = _parse_normalized_matrix(lines)
    assert len(X_norm) == X.nrow()
    if X_norm:
        assert len(X_norm[0]) == X.ncol()

    Y = X.matrix()
    Y._X = X_norm
    return Y

def _format_exec_file(X_file, m_file, binreg_path):
    from StringIO import StringIO

    handle = StringIO()
    w = handle.write
    w("X = load('%s');\n" % X_file)
    w("m = load('%s');\n" % m_file)
    w("m = m(1,:);\n")   # convert m from matrix to vector.
    w("\n");
    w("curpath = cd;\n")
    w("cd('%s');\n" % binreg_path)
    w("Xout = quantnorm(X, m);\n")
    w("cd(curpath);\n")
    w("\n")
    w("disp('NORMALIZED MATRIX')\n")
    w("disp(num2str(Xout, 16))\n")
    w("quit;\n")
    handle.seek(0)
    return handle.read()

