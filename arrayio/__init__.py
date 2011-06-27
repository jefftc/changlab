"""

choose_format   Format to read a file.
guess_format    Format for a matrix.
convert

read
write

"""
# TODO:
# Should detect available FORMATS automatically.

import os

from const import *

FORMAT_NAMES = [
    # Most specific to most general format.
    #"res_format",     # Characteristic row lengths.  Matrix specific.
    "gct_format",     # Matrix more specific than ODF format.
    #"odf_format",     # Very specific header.
    "jeffs_format", 
    "cdt_format", 
    "pcl_format", 
    "tab_delimited_format",
    "csv_format",
    ]
FORMATS = [__import__(x, globals(), locals(), []) for x in FORMAT_NAMES]
tdf = tab_delimited_format  # for convenience

def choose_format(locator):
    # Return the module that can read this format or None.
    for format in FORMATS:
        if format.is_format(locator):
            return format
    return None

def guess_format(X):
    # Return a module that describes a matrix or None.
    for format in FORMATS:
        if format.is_matrix(X):
            return format
    return None

def read(locator, datatype=float, format=None):
    # Bug: this function fails if passed a file handle.
    from genomicode import filelib
    format = format or choose_format(locator)
    if format is None:
        msg = "I could not find a format for: %r" % locator
        if not os.path.exists(locator):
            msg = "%s.  It does not appear to be a file." % msg
        assert format is not None, msg
    #print "PARSING WITH", format
    return format.read(locator, datatype=datatype)

def write(X, handle, format=None):
    format = format or guess_format(X)
    assert format is not None
    format.write(X, handle)

## def _res_to_gct(X):
##     # Will lose the CALL information.
##     # Will lose the SCALE_FACTOR.
##     # Will lose the column DESCRIPTIONS.
##     from genomicode import Matrix
##     assert res_format.is_matrix(X)

##     # Figure out the annotation names for the row name and description.
##     acc_header, desc_header, call_header = X.row_headers()
##     if X._synonyms[ROW_ID] != acc_header:
##         acc_header, desc_header = desc_header, acc_header
##     assert X._synonyms[ROW_ID] == acc_header

##     row_names = X._row_names
##     col_names = X._col_names
##     row_headers = ["Name", "Description"]
##     col_headers = None
##     row_annots = {}
##     col_annots = {}
##     synonyms = { ROW_ID : "Name" }

##     row_annots["Name"] = X.row_annots(acc_header)
##     row_annots["Description"] = X.row_annots(desc_header)

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert gct_format.is_matrix(x)
##     return x
    
## def _res_to_odf(X):
##     # Will lose the CALL information.
##     # Will lose the SCALE_FACTOR.
    
##     # Same as the GCT converter, except for the col_annots.
##     x = _res_to_gct(X)
##     assert not x._col_annots and not x._col_headers
##     assert len(x._row_annots) == 2
##     x._col_annots["TYPE"] = ["float"]*x.ncol()
##     x._col_annots["DESCRIPTION"] = X.col_annots("DESCRIPTION")
##     x._col_headers = sorted(x._col_annots)
##     assert odf_format.is_matrix(x)
##     return x
    
## def _res_to_pcl(X):
##     # Will lose the CALL information.
##     # Will lose the SCALE_FACTOR.
##     # Will lose the column DESCRIPTIONS.
##     from genomicode import Matrix
##     assert res_format.is_matrix(X)
    
##     # Figure out the annotation names for the row name and description.
##     acc_header, desc_header, call_header = X.row_headers()
##     if X._synonyms[ROW_ID] != acc_header:
##         acc_header, desc_header = desc_header, acc_header
##     assert X._synonyms[ROW_ID] == acc_header

##     # Make sure the names don't conflict.
##     row_name = acc_header
##     if row_name == "NAME":
##         row_name = "NAME0"

##     row_names = X._row_names
##     col_names = X._col_names
##     row_annots = {}
##     col_annots = {}
##     row_headers = [row_name, "NAME"]
##     col_headers = None
##     synonyms = { ROW_ID : row_name }
    
##     row_annots[row_name] = X._row_annots[acc_header]
##     row_annots["NAME"] = X._row_annots[desc_header]

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert pcl_format.is_matrix(x)
##     return x
    
## def _gct_to_odf(X):
##     from genomicode import Matrix
##     assert gct_format.is_matrix(X)

##     assert len(X._synonyms) == 1

##     row_names = X._row_names
##     col_names = X._col_names
##     row_headers = X._row_headers
##     col_headers = None
##     row_annots = {}
##     col_annots = {}
##     synonyms = X._synonyms

##     row_annots = X._row_annots
##     col_annots["TYPE"] = ["float"]*X.ncol()

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert odf_format.is_matrix(x)
##     return x
    
def _gct_to_pcl(X):
    from genomicode import Matrix
    
    assert gct_format.is_matrix(X)
    assert len(X.col_names()) == 1
    assert X._col_order and X._col_order[0] == tab_delimited_format.SAMPLE_NAME
    assert len(X.row_names()) == 2
    name, desc = X.row_names()

    row_order = ["GeneID", "NAME"]
    col_order = X._col_order[:]
    row_names = {}
    col_names = X._col_names.copy()
    synonyms = {}

    row_names["GeneID"] = X._row_names[name]
    row_names["NAME"] = X._row_names[desc]
    synonyms[ROW_ID] = "GeneID"
    synonyms[COL_ID] = col_order[0]

    x = Matrix.InMemoryMatrix(
        X._X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    x = Matrix.add_synonyms(x, synonyms)
    assert pcl_format.is_matrix(x)
    return x
    
## def _odf_to_gct(X):
##     # Will lose the column annotations.
##     from genomicode import Matrix
##     assert odf_format.is_matrix(X)

##     # Figure out the annotation names for the row name and description.
##     x = X.row_headers()
##     assert len(x) == 2
##     row_name, row_desc = x
##     if X._synonyms[ROW_ID] != row_name:
##         row_name, row_desc = row_desc, row_name
##     assert X._synonyms[ROW_ID] == row_name

##     row_names = X._row_names
##     col_names = X._col_names
##     row_headers = ["NAME", "DESCRIPTION"]
##     col_headers = None
##     row_annots = {}
##     col_annots = {}
##     synonyms = { ROW_ID : "NAME" }
    
##     row_annots["NAME"] = X._row_annots[row_name]
##     row_annots["DESCRIPTION"] = X._row_annots[row_desc]

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert gct_format.is_matrix(x)
##     return x

## def _odf_to_pcl(X):
##     # Will lose the column annotations.
##     from genomicode import Matrix
##     assert odf_format.is_matrix(X)

##     # Figure out the annotation names for the row name and description.
##     x = X.row_headers()
##     assert len(x) == 2
##     row_name, row_desc = x
##     if X._synonyms[ROW_ID] != row_name:
##         row_name, row_desc = row_desc, row_name

##     # Make sure the names don't conflict.
##     new_row_name = row_name
##     if new_row_name == "NAME":
##         new_row_name = "NAME0"

##     row_names = X._row_names
##     col_names = X._col_names
##     row_headers = [new_row_name, "NAME"]
##     col_headers = None
##     row_annots = {}
##     col_annots = {}
##     synonyms = { ROW_ID : new_row_name }
    
##     row_annots[new_row_name] = X._row_annots[row_name]
##     row_annots["NAME"] = X._row_annots[row_desc]

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert pcl_format.is_matrix(x)
##     return x

def _jeff_to_gct(X):
    # Will lose a lot of annotations.
    from genomicode import Matrix
    assert jeffs_format.is_matrix(X)

    assert len(X.col_names()) == 1
    assert X._col_order and X._col_order[0] == tab_delimited_format.SAMPLE_NAME

    row_order = ["NAME", "DESCRIPTION"]
    col_order = X._col_order[:]
    row_names = {}
    col_names = X._col_names.copy()
    synonyms = {}

    row_names["NAME"] = X._row_names["Probe.Set.ID"]
    row_names["DESCRIPTION"] = X._row_names["Gene.Symbol"]
    synonyms[ROW_ID] = "NAME"
    synonyms[COL_ID] = col_order[0]

    x = Matrix.InMemoryMatrix(
        X._X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    x = Matrix.add_synonyms(x, synonyms)
    assert gct_format.is_matrix(x)
    return x

## def _jeff_to_odf(X):
##     # Will lose a lot of annotations.
##     # Basically like GCT format, but with a TYPE column annotation.
##     x = _jeff_to_gct(X)
    
##     assert not x._col_annots and not x._col_headers
##     assert len(x._row_annots) == 2
##     x._col_annots["TYPE"] = ["float"]*x.ncol()
##     x._col_headers = sorted(x._col_annots)
##     assert odf_format.is_matrix(x)
##     return x

def _jeff_to_pcl(X):
    # Will lose a lot of annotations.
    from genomicode import Matrix
    assert jeffs_format.is_matrix(X)

    assert len(X.col_names()) == 1
    assert X._col_order and X._col_order[0] == tab_delimited_format.SAMPLE_NAME

    row_order = ["Probe.Set.ID", "NAME"]
    col_order = X._col_order[:]
    row_names = {}
    col_names = X._col_names.copy()
    synonyms = {}

    row_names["Probe.Set.ID"] = X._row_names["Probe.Set.ID"]
    row_names["NAME"] = X._row_names["Gene.Symbol"]
    synonyms[ROW_ID] = "Probe.Set.ID"
    synonyms[COL_ID] = col_order[0]

    x = Matrix.InMemoryMatrix(
        X._X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    x = Matrix.add_synonyms(x, synonyms)
    assert pcl_format.is_matrix(x)
    return x

def _pcl_to_gct(X):
    # Will lose the column annotations.
    from genomicode import Matrix
    from genomicode import parselib
    assert pcl_format.is_matrix(X)

    assert len(X.col_names()) == 1
    assert X._col_order and X._col_order[0] == tab_delimited_format.SAMPLE_NAME
    assert len(X.row_names()) > 0

    # Figure out the annotation names for the row name and description.
    if len(X.row_names()) == 1:
        row_name, row_desc = X.row_names()[0], None
    else:
        row_name, row_desc = X.row_names()[:2]
    
    row_order = ["NAME", "DESCRIPTION"]
    col_order = X._col_order[:]
    row_names = {}
    col_names = X._col_names.copy()
    synonyms = {}

    row_names["NAME"] = X._row_names[row_name]
    if row_desc:
        row_names["DESCRIPTION"] = X._row_names[row_desc]
    else:
        # Make up default row names.
        x = ["DESC%s" % x for x in parselib.pretty_range(0, X.nrow())]
        row_names["DESCRIPTION"] = x
    synonyms[ROW_ID] = "NAME"
    synonyms[COL_ID] = col_order[0]

    x = Matrix.InMemoryMatrix(
        X._X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    x = Matrix.add_synonyms(x, synonyms)
    assert gct_format.is_matrix(x)
    return x

## def _pcl_to_odf(X):
##     # Will lose GWEIGHT, GORDER, EWEIGHT, EORDER.
##     from genomicode import Matrix
##     assert pcl_format.is_matrix(X)

##     assert len(X._synonyms) == 1

##     row_names = X._row_names
##     col_names = X._col_names
##     row_headers = X._row_headers[:2]
##     col_headers = None
##     row_annots = {}
##     col_annots = {}
##     synonyms = X._synonyms

##     row_annots = X._row_annots
##     col_annots["TYPE"] = ["float"]*X.ncol()

##     x = Matrix.InMemoryMatrix(
##         X._X, row_names=row_names, col_names=col_names,
##         row_headers=row_headers, col_headers=col_headers,
##         row_annots=row_annots, col_annots=col_annots, synonyms=synonyms)
##     assert odf_format.is_matrix(x)
##     return x
    
def _tdf_to_gct(X):
    from genomicode import Matrix
    from genomicode import parselib
    assert tab_delimited_format.is_matrix(X)

    assert len(X.col_names()) == 1
    assert X._col_order and X._col_order[0] == tab_delimited_format.SAMPLE_NAME

    name_header = "NAME"
    desc_header = "DESCRIPTION"

    # Make up default names.
    name = ["NAME%s" % x for x in parselib.pretty_range(0, X.nrow())]
    desc = ["DESC%s" % x for x in parselib.pretty_range(0, X.nrow())]
    
    # Try to find better names.
    if not X.row_names():
        pass
    elif len(X.row_names()) == 1:
        # Only 1 header, so use that for the name.
        name = X.row_names(X.row_names()[0])
    else:
        # Use the first two columns for the name and description.
        name_i, desc_i = 0, 1
        # See if there is a ROW_ID set.  If there is, use that for NAME.
        if hasattr(X, "_synonyms") and ROW_ID in X._synonyms:
            name_i = X.row_names().index(X._synonyms[ROW_ID])
            if name_i == desc_i:
                # name_i used to be 0, and desc_i is not 0.
                assert desc_i != 0
                desc_i = 0
        assert name_i != desc_i
        name = X.row_names(X.row_names()[name_i])
        desc = X.row_names(X.row_names()[desc_i])
    
    row_order = [name_header, desc_header]
    col_order = X._col_order[:]
    row_names = {}
    col_names = X._col_names.copy()
    synonyms = {}
    
    row_names[name_header] = name
    row_names[desc_header] = desc
    synonyms[ROW_ID] = name_header
    synonyms[COL_ID] = col_order[0]
        
    x = Matrix.InMemoryMatrix(
        X._X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    x = Matrix.add_synonyms(x, synonyms)
    assert gct_format.is_matrix(x)
    return x

def _any_to_csv(X):
    return X

def _any_to_tdf(X):
    return X

def _format2name(format):
    name = format.__name__
    if "." in name:
        name = name[name.rindex(".")+1:]
    return name

def _choose_converter(from_format, to_format):
    # Get the names of the formats.
    from_name = _format2name(from_format)
    to_name = _format2name(to_format)

    # o Doesn't make sense to convert to CDT format.  Must be produced
    #   by a clustering algorithm.
    # o Doesn't make sense to convert to RES format.  Must be produced
    #   by a algorithm that generates CALLs.

    # See if a converter exists.
##     if from_name == "res_format" and to_name == "gct_format":
##         return _res_to_gct
##     elif from_name == "res_format" and to_name == "odf_format":
##         return _res_to_odf
##     elif from_name == "res_format" and to_name == "pcl_format":
##         return _res_to_pcl
##     elif from_name == "gct_format" and to_name == "odf_format":
##         return _gct_to_odf
    if from_name == "gct_format" and to_name == "pcl_format":
        return _gct_to_pcl
##     elif from_name == "odf_format" and to_name == "gct_format":
##         return _odf_to_gct
##     elif from_name == "odf_format" and to_name == "pcl_format":
##         return _odf_to_pcl
    elif from_name == "jeffs_format" and to_name =="gct_format":
        return _jeff_to_gct
##     elif from_name == "jeffs_format" and to_name == "odf_format":
##         return _jeff_to_odf
    elif from_name == "jeffs_format" and to_name == "pcl_format":
        return _jeff_to_pcl
    elif from_name == "pcl_format" and to_name == "gct_format":
        return _pcl_to_gct
##     elif from_name == "pcl_format" and to_name == "odf_format":
##         return _pcl_to_odf
    elif from_name == "tab_delimited_format" and to_name == "gct_format":
        return _tdf_to_gct
    elif to_name == "csv_format":
        return _any_to_csv
    elif to_name == "tab_delimited_format":
        return _any_to_tdf

    # No converter exists.
    return None

def convert(X, from_format=None, to_format=None):
    if from_format is None:
        from_format = guess_format(X)
        assert from_format is not None, "Could not find format for matrix X."
    if to_format is None:
        to_format = "tab_delimited_format"
    if from_format is to_format:
        return X
    convert_fn = _choose_converter(from_format, to_format)
    assert convert_fn is not None, \
           "I could not find a converter from %s to %s." % (
        _format2name(from_format), _format2name(to_format))
    return convert_fn(X)
