"""

Functions:
normalize

"""
import os, sys

def normalize(X, Y, matlab=None, binreg_path=None):
    # X and Y are Matrixes.  Will normalize X to Y and return a
    # normalized version of X.
    #import copy
    import subprocess
    import tempfile
    import binreg
    from dwdnorm import _write_matlab_matrix, _parse_normalized_matrix
    from dwdnorm import _safe_unlink
    
    # Set defaults.
    matlab = matlab or "matlab"
    binreg_path = binreg.find_binreg_20(binreg_path)
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
    X_file = Y_file = None
    try:
        x, X_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        x, Y_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        _write_matlab_matrix(X, X_file)
        _write_matlab_matrix(Y, Y_file)
        script = _format_exec_file(X_file, Y_file, binreg_path)
        w.write(script)
        w.close()
        lines = r.readlines()
    finally:
        _safe_unlink(X_file)
        _safe_unlink(Y_file)

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

def _format_exec_file(X_file, Y_file, binreg_path):
    from StringIO import StringIO
    #import dwdnorm

    #X_mat = dwdnorm._format_matlab_matrix(X)
    #Y_mat = dwdnorm._format_matlab_matrix(Y)
    handle = StringIO()
    w = handle.write
    #w("X = %s;\n" % X_mat)
    #w("Y = %s;\n" % Y_mat)
    w("X = load('%s');\n" % X_file)
    w("Y = load('%s');\n" % Y_file)
    w("\n");
    w("curpath = cd;\n")
    w("cd('%s');\n" % binreg_path)
    w("Xout = std_rows_to_y(X, Y);\n")
    w("cd(curpath);\n")
    w("\n")
    w("disp('NORMALIZED MATRIX')\n")
    w("disp(num2str(Xout, 16))\n")
    #w("save('%s', 'Xout', '-ASCII', '-TABS');\n" % outfile)
    w("quit;\n")
    handle.seek(0)
    return handle.read()

