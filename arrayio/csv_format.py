"""Generic comma-delimited format.

Same format as tab_delimited_format, but uses commas as delimiter
instead of tabs.


Functions:
read
write
is_format
is_matrix

"""
import os, sys

def is_format(locator_str):
    from genomicode import filelib
    if not filelib.exists(locator_str):
        return False
    handle = filelib.openfh(locator_str)
    x = handle.readline()
    handle.close()   # need to close it properly, or gunzip might not die.
    if not x:  # blank file
        return False
    if "," in x:
        return True
    return False

def is_matrix(X):
    # Any matrix can be csv format.
    return True

def read(handle, hrows=None, hcols=None, datatype=float):
    import csv
    from StringIO import StringIO
    from genomicode import filelib
    import tab_delimited_format

    # Convert this to tab-delimited format and let the other module
    # deal with it.
    outhandle = StringIO()
    reader = csv.reader(filelib.openfh(handle))
    for row in reader:
        print >>outhandle, "\t".join(row)
    outhandle.seek(0)
    return tab_delimited_format.read(
        outhandle, hrows=hrows, hcols=hcols, datatype=datatype)

def write(X, handle):
    # Who wants to write CSV?
    # If implement, need to be careful to escape commas in gene
    # descriptions, etc.
    raise NotImplementedError
