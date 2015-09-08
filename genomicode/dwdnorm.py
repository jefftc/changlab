"""

Functions:
find_dwd
normalize

"""
# _format_matlab_matrix    XXX should move to matlab.py library.
# _format_matlab_vector
# _write_matlab_matrix

# _format_exec_file
# _parse_normalized_matrix

import os, sys

def find_dwd(default_path):
    # Return the path to BatchAdjust or None.
    import config
    
    search_paths = [
        default_path,
        config.BatchAdjust,
        "BatchAdjust",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        # Test for some well-known BatchAdjust files.
        files = [
            "ReadMe.txt", "SubRoutines", "basm1.m",
            # Files removed from distribution to save space.
            #"Data", "BatchAdjustSM.m", "BatchAdjustSMtest.m", 
            ]
        complete = True   # Is this distribution complete.
        for file in files:
            filename = os.path.join(spath, file)
            if not os.path.exists(filename):
                complete = False
                break
        if not complete:
            continue
        path = spath
        break
    return path


def normalize(X, Y, version=None, matlab=None, dwd_path=None):
    # X is a Matrix of the data.  Y is a vector of -1, 1 indicating
    # the class of each sample.  Returns the normalized version of the
    # X matrix.  version is "default" or "bild", which indicates how
    # to handle the normalization.  If version is "bild", then Y can
    # consist of -1, 1, or 2, where -1 and 1 distinguish the two
    # classes of the training set, and 2 indicates the validation set.
    # Returns the normalized Matrix.
    #import copy
    import subprocess
    import tempfile
    import Matrix
    
    # Set defaults.
    version = version or "default"
    matlab = matlab or "matlab"
    dwd_path = find_dwd(dwd_path)
    assert dwd_path, "I could not find DWD."
    dwd_path = os.path.realpath(dwd_path)
    temp_path = "."
    assert version in ["default", "bild"]

    # Start and instance of matlab.
    matlab_args = [
        "-nosplash", "-nodesktop", "-nodisplay", "-nojvm"]
    x = " ".join(matlab_args)
    cmd = "%s %s" % (matlab, x)
    #w, r = os.popen4(cmd, bufsize=0)
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout

    # Run the normalization.
    X_file = None
    try:
        x, X_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        _write_matlab_matrix(X, X_file)
        script = _format_exec_file(X_file, Y, version, dwd_path)
        # This might fail with a "Broken pipe" if no Matlab installed.
        w.write(script)
        w.close()
        lines = r.readlines()
    finally:
        _safe_unlink(X_file)

    #open("dwd_norm.m", 'w').write(script)
    X_norm = _parse_normalized_matrix(lines)
    assert len(X_norm) == X.nrow()
    if X_norm:
        assert len(X_norm[0]) == X.ncol()

    #Y = Matrix.InMemoryMatrix(
    #    X_norm, row_names=X._row_names, col_names=X._col_names,
    #    row_headers=X._row_headers, col_headers=X._col_headers,
    #    row_annots=X._row_annots, col_annots=X._col_annots,
    #    synonyms=X._synonyms)
    #Y = copy.deepcopy(X)
    Y = X.matrix()
    Y._X = X_norm
    return Y
    
def _format_matlab_matrix(X):
    X = X.slice()
    lines = [" ".join(map(str, x)) for x in X]
    return "[%s]" % "\n".join(lines)

def _format_matlab_vector(X):
    return "[%s]" % " ".join(map(str, X))

def _write_matlab_matrix(X, filename):
    import Matrix
    assert Matrix.is_Matrix(X), "X must be a Matrix object."
    
    handle = open(filename, 'w')
    X = X.slice()
    for x in X:
        print >>handle, " ".join(map(str, x))
    handle.close()

def _format_exec_file(X_file, Y, version, dwd_path):
    from StringIO import StringIO

    assert version in ["default", "bild"]
    function = "BatchAdjustSM"
    if version == "bild":
        function = "basm1"

    #X_mat = _format_matlab_matrix(X)
    Y_mat = _format_matlab_vector(Y)

    handle = StringIO()
    w = handle.write
    w("params = struct( ...\n")
    w("  'viplot', zeros(4, 1), ...\n")
    w("  'savestr', 'out.dat', ...\n")
    w("  'iscreenwrite', 1);\n")
    w("\n")
    w("X = load('%s');\n" % X_file)
    w("Y = %s;\n" % Y_mat)
    w("\n");
    w("curpath = cd;\n")
    w("cd('%s');\n" % dwd_path)
    w("Xout = %s(X, Y, params);\n" % function)
    w("cd(curpath);\n")
    w("\n")
    w("disp('NORMALIZED MATRIX')\n")
    w("disp(num2str(Xout, 16))\n")
    #w("save('%s', 'Xout', '-ASCII', '-TABS');\n" % outfile)
    w("quit;\n")
    handle.seek(0)
    x = handle.read()
    return x

def _safe_unlink(filename):
    if filename and os.path.exists(filename):
        os.unlink(filename)

def _parse_normalized_matrix(lines):
    # Read and parse the results.
    # Format:
    # <BatchAdjust output>
    # >> >> >> >> >> >> >> >> >> NORMALIZED MATRIX
    # >> <row 0>
    # <row 1>
    # ...
    # >>
    # >> Warning: failed to create preference directory /home/jchang/.matlab/R2
    # Check directory permissions.
    # 
    # Pull out every line after NORMALIZED MATRIX.
    # May have a warning at the end.  Ignore it.
    for i in range(len(lines)):
        if lines[i].find("NORMALIZED MATRIX") >= 0:
            break
    else:
        raise AssertionError, "I could not find the normalized output."
    lines = lines[i+1:]

    # Parse out the matrix.
    X_norm = []
    for line in lines:
        line = line.replace(">", "").strip()
        if not line:
            continue
        if line.startswith("Warning: failed to create preference directory"):
            continue
        if line.startswith("Check directory permissions."):
            continue
        x = [float(x) for x in line.split()]
        #assert len(x) == X.ncol(), "%d %d" % (len(x), X.ncol())
        X_norm.append(x)
        assert len(X_norm[0]) == len(x)
    return X_norm
    
