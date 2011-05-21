"""

Functions:
split_tdf   Split a tab-delimited format
read_tdf
write_tdf

strip_each

safe_float  Moved to jmath.
safe_int    Moved to jmath.

"""
import os, sys

def split_tdf(data):
    # Bug: will skip blank lines.
    return [x.rstrip("\r\n").split("\t") for x in data.split("\n")
            if x.strip()]

def read_tdf(handle, num_header_rows, num_header_cols):
    # Format:
    # - First num_header_cols contain metadata, last ones samples.
    # - First row headings, last ones samples.
    # - Return as a single matrix.
    import filefns
    import jmath

    handle = filefns.openfh(handle)
    data = split_tdf(handle.read())

    X = []
    num_cols = None
    for cols in data:
        if num_cols is None:
            num_cols = len(cols)
        assert len(cols) == num_cols
        if len(X) >= num_header_rows:
            x = map(jmath.safe_float, cols[num_header_cols:])
            cols[num_header_cols:] = x
        X.append(cols)
    return X

def write_tdf(data, outhandle):
    for x in data:
        x = map(str, x)
        print >>outhandle, "\t".join(x)

def strip_each(L):
    return [x.strip() for x in L]


# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    import ciofns
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in ciofns.__dict__.keys():
        if not name.startswith("__"):
            this_module.__dict__[name] = ciofns.__dict__[name]
