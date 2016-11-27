"""

Functions:
split_tdf   Split a tab-delimited format
read_tdf    OBSOLETE
write_tdf   OBSOLETE

strip_each

"""
import os, sys

def split_tdf(data, strip=False):
    # Bug: will skip blank lines.
    matrix = [
        x.rstrip("\r\n").split("\t") for x in data.split("\n") if x.strip()]
    if strip:
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = matrix[i][j].strip()
    return matrix

## def read_tdf(handle, num_header_rows, num_header_cols):
##     # Format:
##     # - First num_header_cols contain metadata, last ones samples.
##     # - First row headings, last ones samples.
##     # - Return as a single matrix.
##     import filelib
##     import jmath

##     handle = filelib.openfh(handle)
##     data = split_tdf(handle.read())

##     X = []
##     num_cols = None
##     for cols in data:
##         if num_cols is None:
##             num_cols = len(cols)
##         assert len(cols) == num_cols
##         if len(X) >= num_header_rows:
##             x = map(jmath.safe_float, cols[num_header_cols:])
##             cols[num_header_cols:] = x
##         X.append(cols)
##     return X

## def write_tdf(data, outhandle):
##     for x in data:
##         x = map(str, x)
##         print >>outhandle, "\t".join(x)

CLEAN_RE = None
CLEAN_DISALLOWED = None
def _py_cleanwrite(data, outhandle, delim="\t"):
    global CLEAN_RE
    global CLEAN_DISALLOWED
    import re

    disallowed = "\r\n" + delim
    if CLEAN_RE is None or CLEAN_DISALLOWED != disallowed:
        CLEAN_RE = re.compile("[%s]" % disallowed)
        CLEAN_DISALLOWED = disallowed

    for x in data:
        x = x[:]
        for i in range(len(x)):
            if x[i] is None:
                x[i] = ""
        x = map(str, x)
        x = [CLEAN_RE.subn(" ", x)[0].strip() for x in x]
        x = delim.join(x)
        print >>outhandle, x
_cleanwrite = _py_cleanwrite

def cleanwrite(data, outhandle, delim="\t"):
    import types
    
    # The C version of _cleanwrite requires outhandle to be a file
    # object and won't work with StringIO objects.  If the user does
    # not provide a real file object, then use the python
    # implementation.
    if type(outhandle) is not types.FileType:
        _py_cleanwrite(data, outhandle, delim=delim)
        return
    #print "Calling _cleanwrite"
    _cleanwrite(data, outhandle, delim=delim)

def strip_each(L):
    return [x.strip() for x in L]


# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    #raise ImportError
    import ciolib
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in ciolib.__dict__.keys():
        if not name.startswith("__"):
            this_module.__dict__[name] = ciolib.__dict__[name]
