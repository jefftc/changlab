"""

Stores a 2D matrix of objects in the format of:
                <ROW_HEADER_0> ...  <name_0>  ...  <name_n>  
<COL_HEADER_0>                      <name_0>  ...  <name_n>
     ...                              ...     ...    ...
  <name_0>         <name_0>    ...    X_00    ...    X_0N
  <name_1>         <name_1>    ...    X_10    ...    X_1N
    ...              ...       ...    ...     ...    ...
  <name_n>         <name_n>    ...    X_N0    ...    X_NN

There is a matrix of mxn values X where the rows and columns have
optional headers.  The first row and column contains headers that have
no names.  How to handle this is up to the user of the Matrix.


X = Matrix(...)

X.dim()                num rows, num cols
X.nrow()
X.ncol()
X.row_names()          Return the headers.
X.row_names(header)    Return the names under header.


X.matrix(...) -> <Matrix>               Parameters same as X.value.
X.slice(...) -> <list> of <list>s       Parameters same as X.value.
X.value(<int>, <int>) ->  <obj>         <int> indexes can be negative.
X.value(<int>, None) -> <list>
X.value(<int>, <list>) -> <list>
X.value(<int>, (<int>, <int>)) -> <list>
X.value(<int>, (<int>, None)) -> <list>
X.value(<list>, <list>) -> <list> of <list>s
X.value(row=<str>, row_header=<str>) -> <list>
X.value(row=<str>, row_header=<list>) -> <list> of <list>s
X.value(row=<str>, row_header=None) -> <list> of <list>s
X.value(row=<list>, row_header=<str>) -> <list>         #  <list> of <list>s???
X.value(row=<list>, row_header=<list>) -> <list> of <list>s
X.value(row=<list>, row_header=None) -> <list> of <list>s
X.value(row=<list>, row_header=<list>, col=<list>, col_annot=<list>) ->
  list> of <list>s

X[(<int>, <int>)]                       Shorthand for X.value
X[(<int>, None)]
X[(<int>, (<int>, <int>))]
X[(<int>, <list>)]
X[(<list>, <list>)]


Classes:
AbstractMatrix
InMemoryMatrix

"""
# TODO
# - allow case insensitive matching of annotations

import os, sys

class AbstractMatrix:
    """

    Methods:
    matrix
    slice
    value

    dim
    ncol
    nrow
    
    row_names
    col_names
    
    """
    def __init__(self, synonyms):
        # Set some variables for caching.
        self._row_header2name2indexes = {}  # header -> name -> list of indexes
        self._col_header2name2indexes = {}  # header -> name -> list of indexes
        assert type(synonyms) is type({})
        self._synonyms = synonyms

    def _resolve_synonym(self, header, name_fn, synonyms):
        from types import StringType, TupleType, ListType

        # Sets introduced in Python 2.4.
        if sys.hexversion >= 0x020400F0:
            if type(header) is type(set()):
                header = list(header)

        # header can be None, string, or list of strings.
        if type(header) is StringType:
            if header not in name_fn() and header in synonyms:
                header = synonyms[header]
        elif type(header) in [TupleType, ListType]:
            for i in range(len(header)):
                if header[i] not in name_fn() and header[i] in synonyms:
                    header[i] = synonyms[header[i]]
        else:
            assert header is None
        return header

    def _index(self, row=None, col=None, row_header=None, col_header=None):
        # Return a tuple (I_row, I_col) where I_row are the indexes of
        # the rows to take, and I_col are the indexes of the cols to
        # take.  I_row and I_col can be None, which indicates taking
        # all indexes.
        row_header = self._resolve_synonym(
            row_header, self.row_names, self._synonyms)
        col_header = self._resolve_synonym(
            col_header, self.col_names, self._synonyms)

        #if not hasattr(self, "_count"):
        #    self._count = 0
        #assert self._count != 4
        #self._count += 1
        
        # This method cannot be used for retrieving annotations,
        # without specifying any ids.
        assert not (row_header and row is None), "must specify row ids"
        assert not (col_header and col is None), "must specify col ids"

        I_row = self._index_one(
            row_header, row, 
            self.row_names, self._get_row_header2name2indexes, self.nrow())
        # Optimization.  Most of the time, the user wants to index columns.
        if col_header is None and col is None:
            I_col = None
        else:
            I_col = self._index_one(
                col_header, col, 
                self.col_names, self._get_col_header2name2indexes, self.ncol())
        return I_row, I_col

    def _index_one(
        self, header, name, get_headers_fn, get_name2indexes_fn, length):
        # Return a list of indexes.
        import operator
        from types import IntType, StringType, TupleType, ListType

        # Cases to handle:
        #     name             header  Notes
        #  1  None             None    Return None (all indexes).
        #  2  <int>            None    <int> can be negative.
        #  3  <list> of <int>  None    Each of the indexes can be negative.
        #  4  (<s>, <e>)       None    <s> and <e> can be negative.
        #  5  <str>            None    Indexes that match any header.
        #  6  <str>            <str>   Indexes that match a specific header.
        #  7  <str>            <list>  Indexes that match the listed headers.
        #  8  <list> of <str>  None    Indexes that match any header.
        #  9  <list> of <str>  <str>   Indexes that match a specific headers.
        # 10  <list> of <str>  <list>  Indexes that match the listed headers.
        # 11  empty <list>             Empty (no indexes).

        # Sets introduced in Python 2.4.
        if sys.hexversion >= 0x020400F0:
            if type(name) is type(set()):
                name = list(name)
            if type(header) is type(set()):
                header = list(header)
        type_name = type(name)

        case = None
        # Case 1.
        # For slicing, making copies.
        if name is None and header is None:
            case = 1
            indexes = None
        # Case 9.
        # Used very often for selecting genes.
        elif (type_name is ListType and list_type(name) is StringType
              and type(header) is StringType):
            case = 9
            indexes = self._find_names(name, [header], get_name2indexes_fn)
        # Case 2.
        elif type_name is IntType:
            assert header is None
            case = 2
            indexes = [normalize_index(name, length)]
        # Case 3.
        elif (type_name is ListType and list_type(name) is type(0)
              and header is None):
            case = 3
            indexes = normalize_indexes(name, length)
        # Case 4.
        elif type_name is TupleType:
            assert len(name) == 2, "Must supply a start and end index."
            assert header is None
            case = 4
            indexes = normalize_slice(name[0], name[1], length)
        # Case 5.        
        elif type_name is StringType and header is None:
            case = 5
            headers = get_headers_fn()
            indexes = self._find_names([name], headers, get_name2indexes_fn)
        # Case 6.
        elif type_name is StringType and type(header) is StringType:
            case = 6
            indexes = self._find_names([name], [header], get_name2indexes_fn)
        # Case 7.
        elif type_name is StringType and type(header) is ListType:
            case = 7
            indexes = self._find_names([name], header, get_name2indexes_fn)
        # Case 8.
        elif (type_name is ListType and list_type(name) is StringType
              and header is None):
            case = 8
            headers = get_headers_fn()
            indexes = self._find_names(name, headers, get_name2indexes_fn)
        # Case 10.
        elif (type_name is ListType and list_type(name) is StringType
              and type(header) is ListType):
            case = 10
            indexes = self._find_names(name, header, get_name2indexes_fn)
        # Case 11.
        elif (type_name is ListType and not len(name)):
            case = 11
            indexes = []
        else:
            raise AssertionError, "I don't know how to interpret the indexes."
        #print "CASE %d" % case   # DEBUG

        return indexes

    def _get_header2name2indexes(self, header, header2name2indexes, names_fn):
        if header not in header2name2indexes:
            # This line is only 5% of the running time.
            names = names_fn(header)
            x = _index_strings(names)
            header2name2indexes[header] = x
        return header2name2indexes[header]

    def _get_row_header2name2indexes(self, header):
        return self._get_header2name2indexes(
            header, self._row_header2name2indexes, self.row_names)

    def _get_col_header2name2indexes(self, header):
        return self._get_header2name2indexes(
            header, self._col_header2name2indexes, self.col_names)

    def _find_names(self, names, headers, get_name2indexes_fn):
        # Look for the names in the appropriate headers and return
        # their indexes.  The indexes will be in the same order as the
        # names, but they will not necessarily be parallel due to
        # missing or duplicated names.
        #
        # Currently, a blank name will match anything.  Is this really
        # the desired behavior?
        header2name2indexes = [
            get_name2indexes_fn(header) for header in headers]

        indexes = [None] * max(self.nrow(), self.ncol())
        num = 0
        seen = {}
        for name in names:
            for i, header in enumerate(headers):
                name2indexes = header2name2indexes[i]
                if name not in name2indexes:
                    continue
                # Add the indexes, being careful not to include
                # duplicate ones.
                for i in name2indexes[name]:
                    if i in seen:
                        continue
                    seen[i] = 1
                    indexes[num] = i
                    num += 1
        return indexes[:num]

    def matrix(self, row=None, col=None, row_header=None, col_header=None):
        x = self._index(row, col, row_header, col_header)
        I_row, I_col = x
        return self._matrix(I_row, I_col)

    def slice(self, row=None, col=None, row_header=None, col_header=None):
        x = self._index(row, col, row_header, col_header)
        I_row, I_col = x
        return self._slice(I_row, I_col)

    def value(self, row=None, col=None, row_header=None, col_header=None):
        X = self.slice(row, col, row_header, col_header)

        # Figure out the type of the return value.
        # Cases to handle:
        #  1  Both row and col are <int>s.                 <obj>
        #  4  row is <int>.                                <list>, by row
        #  5  col is <int>.                                <list>, by col
        # 11  Everything else.                             <list> of <list>s
        #
        # BROKEN:
        #  6  row_header is <str> and col_header is None.  <list>, by col
        #  7  col_header is <str> and row_header is None.  <list>, by row
        # 
        # OBSOLETE:
        #  2  row is <int> and col_name is <str>.          <obj>
        #  3  col is <int> and row_name is <str>.          <obj>
        #  8  Both row_name and col_name are <str>s.       <obj>
        #  9  row_name is <str>.                           <list>, by row
        # 10  col_name is <str>.                           <list>, by col

        OBJ, LIST_BY_ROW, LIST_BY_COL, LIST_OF_LISTS = range(4)
        return_type = None
        # Case 1.
        if type(row) is type(0) and type(col) is type(0):
            return_type = OBJ
        # Case 2.
        #elif type(row) is type(0) and type(col_name) is type(""):
        #    return_type = OBJ
        # Case 3.
        #elif type(col) is type(0) and type(row_name) is type(""):
        #    return_type = OBJ
        # Case 4.
        elif type(row) is type(0):
            return_type = LIST_BY_ROW
        # Case 5.
        elif type(col) is type(0):
            return_type = LIST_BY_COL
        # What are cases 6 and 7 for?  I don't know what I was
        # initially thinking.  If row_header is <string>, row is
        # <list>, and col_header is None, then I should get the rows
        # where annot matches row, and all columns.  Shouldn't this be
        # a LIST of LISTs?
        # Case 6.
        #elif type(row_header) is type("") and col_header is None:
        #    return_type = LIST_BY_COL
        # Case 7.
        #elif type(col_header) is type("") and row_header is None:
        #    return_type = LIST_BY_ROW
        # Case 8.
        #elif type(row_name) is type("") and type(col_name) is type(""):
        #    return_type = OBJ
        # Case 9.
        #elif type(row_name) is type(""):
        #    return_type = LIST_BY_ROW
        # Case 10.
        #elif type(col_name) is type(""):
        #    return_type = LIST_BY_COL
        # Case 11.
        else:
            return_type = LIST_OF_LISTS
        assert return_type is not None   # this should not ever happen

        if return_type == OBJ:
            assert len(X) == 1
            assert len(X[0]) == 1
            X = X[0][0]
        elif return_type == LIST_BY_ROW:
            assert len(X) == 1
            X = X[0]
        elif return_type == LIST_BY_COL:
            for x in X:
                assert len(x) == 1
            X = [x[0] for x in X]
        return X

    def __getitem__(self, key):
        # Interpret key.
        # Cases:
        # 1  <int>                                 row
        # 2  <str>                                 Guess name or header.
        # 4  <obj>, <obj>                          row, col
        #
        # OBSOLETE:
        # 3  [<list> of] <str>, [<list> of] <str>  row_name, col_name
        is_two_tuple = type(key) is type(()) and len(key) == 2
        
        list_and_set_types = [type([])]
        # Sets introduced in Python 2.4.
        if sys.hexversion >= 0x020400F0:
            list_and_set_types.append(type(set()))

        params = {}

        # Case 1.
        if type(key) is type(0):
            params["row"] = key
        # Case 2.
        elif type(key) is type(""):
            #if key in self.row_names():
            #    params["row_name"] = key
            #elif key in self.col_names():
            #    params["col_name"] = key
            #if key in self.row_headers():
            if key in self.row_names():
                params["row_header"] = key
            #elif key in self.col_headers():
            elif key in self.col_names():
                params["col_header"] = key
            else:
                raise KeyError, "Unknown name: %s" % key
        # Case 3.
        #elif (is_two_tuple and
        #      (type(key[0]) is type("") or (
        #    type(key[0]) in list_and_set_types and
        #    list_type(key[0]) is type(""))) and
        #      (type(key[1]) is type("") or (
        #    type(key[1]) in list_and_set_types and
        #    list_type(key[1]) is type("")))):
        #    params["row_name"] = key[0]
        #    params["col_name"] = key[1]
        # Case 4.
        elif is_two_tuple:
            params["row"] = key[0]
            params["col"] = key[1]
        else:
            raise AssertionError, "I don't know how to interpret args"
        return self.value(**params)

    # Implement remaining methods in a derived class.
    def _matrix(self, row_indexes, col_indexes):
        # Return another instance of this object.  indexes can be
        # None, which means take the entire row/col.
        raise NotImplementedError
    def _slice(self, row_indexes, col_indexes):
        # Return just the underlying matrix as a list of lists.
        raise NotImplementedError
    
    def dim(self):
        raise NotImplementedError
    def nrow(self):
        return self.dim()[0]
    def ncol(self):
        return self.dim()[1]

    def row_names(self, header=None):
        # Return a list of names.  Should raise a KeyError if
        # the header does not exist.
        raise NotImplementedError, "Please implemented in a derived class."

    def col_names(self, header=None):
        # Return a list of names.  Should raise a KeyError if
        # the header does not exist.
        raise NotImplementedError, "Please implemented in a derived class."
    
    def __str__(self):
        X = self.slice()
        return str(X)
    
class InMemoryMatrix(AbstractMatrix):
    # Members:
    # _X
    # _num_rows
    # _num_cols
    # _row_names
    # _col_names
    # _row_order
    # _col_order
    def __init__(
        self, X, row_names={}, col_names={}, row_order=None, col_order=None,
        synonyms={}):
        # row_names    Dict of header -> list of names.
        # col_names    Dict of header -> list of names.
        # row_order    List of the headers, in order (optional).
        # col_order    List of the headers, in order (optional).
        import math
        
        AbstractMatrix.__init__(self, synonyms)

        # Figure out the size of the matrix.
        nrow = ncol = None
        
        # If X is specified, then figure the nrow and ncol from there.
        # If X is not specified, then either there are no rows or
        # there are no columns.  Try to figure out from the
        # annotations.  If those aren't specified, then there are no
        # rows or columns.
        if X:
            nrow, ncol = len(X), len(X[0])
        if nrow is None and row_names:
            nrow = len(row_names[row_names.keys()[0]])
        if ncol is None and col_names:
            ncol = len(col_names[col_names.keys()[0]])
        nrow, ncol = nrow or 0, ncol or 0

        assert len(X) == nrow, "%d %d" % (len(X), nrow)
        if X:
            assert len(X[0]) == ncol

##         # Generate some default row and column names if necessary.
##         if row_names is None:
##             ndigits = 0
##             if nrow < 10:
##                 # handles nrow==1 border case
##                 ndigits = 1
##             elif nrow >= 10:
##                 ndigits = int(math.ceil(math.log(nrow-1, 10)))
##             row_names = ["R%0*d" % (ndigits, i) for i in range(nrow)]
##         if col_names is None:
##             ndigits = 0
##             if ncol:
##                 ndigits = int(math.ceil(math.log(ncol-1, 10)))
##             col_names = ["C%0*d" % (ndigits, i) for i in range(ncol)]

        # Generate a default order if necessary.
        if row_order is None:
            row_order = sorted(row_names)
        if col_order is None:
            col_order = sorted(col_names)
            
        # Do some basic checking on the input data.
        assert len(X) == nrow
        for i in range(len(X)):
            assert len(X[i]) == ncol
        #assert len(row_names) == nrow,"mismatch %d %d" % (len(row_names),nrow)
        #assert len(col_names) == ncol,"mismatch %d %d" % (len(col_names),ncol)
        assert len(row_names) == len(row_order)
        assert len(col_names) == len(col_order)
        for header, names in row_names.iteritems():
            assert len(names)==nrow, "%s %d %d" % (header, len(names), nrow)
        for header, names in col_names.iteritems():
             assert len(names)==ncol, "%s %d %d" % (header, len(names), ncol)
        assert sorted(row_order) == sorted(row_names)
        assert sorted(col_order) == sorted(col_names)

        for synonym, name in synonyms.iteritems():
            assert name in row_names or name in col_names, "unknown: %s" % name

        # Set the member variables of this object.
        self._X = X
        self._num_rows, self._num_cols = nrow, ncol
        self._row_names = row_names.copy()
        self._col_names = col_names.copy()
        self._row_order = row_order[:]
        self._col_order = col_order[:]

    def row_names(self, header=None):
        if header is None:
            return self._row_order
        header = self._resolve_synonym(header, self.row_names, self._synonyms)
        if header not in self._row_names:
            x = sorted(self._row_names)
            x = ", ".join(x)
            raise KeyError, "%s: headers %s" % (header, x)
        return self._row_names[header]
    def col_names(self, header=None):
        if header is None:
            return self._col_order
        header = self._resolve_synonym(header, self.col_names, self._synonyms)
        return self._col_names[header]
    
    def _matrix(self, row_indexes, col_indexes):
        # Get a subset of the matrix.
        X = self._slice(row_indexes, col_indexes)

        # Get a subset of the annotations.
        rnames = self._row_names
        cnames = self._col_names
        if row_indexes is not None:
            rnames = rnames.copy()
            for header, names in self._row_names.iteritems():
                rnames[header] = [names[i] for i in row_indexes]
        if col_indexes is not None:
            cnames = cnames.copy()
            for header, names in self._col_names.iteritems():
                cnames[header] = [names[i] for i in col_indexes]

        # Special case: for empty matrix, all the rows and column
        # annotations go away.  i.e. if no columns, then this results
        # in an empty matrix, and no need for either row or column
        # annotations.
        if not(X) or not X[0]:
            rnames = rnames.copy()
            cnames = cnames.copy()
            for name in rnames:
                rnames[name] = []
            for name in cnames:
                cnames[name] = []
                
        x = InMemoryMatrix(
            X, row_names=rnames, col_names=cnames, row_order=self._row_order,
            col_order=self._col_order, synonyms=self._synonyms)
        return x
    
    def _slice(self, row_indexes, col_indexes):
        # Return just the underlying matrix.
        # If row_indexes or col_indexes is None, then don't slice the
        # rows or columns.
        # Make a copy of the _X variable.
        X = [x[:] for x in self._X]
        if row_indexes is not None:
            X = [X[i] for i in row_indexes]
        if col_indexes is not None:
            X = [[x[i] for i in col_indexes] for x in X]
        # Don't allow empty rows, e.g. [[]], [[], []].
        X = [x for x in X if x]
        return X

    def dim(self):
        return self._num_rows, self._num_cols


def normalize_index(i, length):
    # Normalize index i so that it ranges from [0, length).
    assert type(i) is type(0), "Invalid index: %r" % i

    # Convert negative index into positive one.
    if i < 0:
        i += length
    if i < 0 or i >= length:
        raise IndexError, "matrix index out of range [%d:%d]" % (i, length)
    return i

def normalize_indexes(L, length):
    return [normalize_index(i, length) for i in L]

def _normalize_slice_h(i, length):
    assert type(i) is type(0)

    if i < -length:
        i = 0
    elif i < 0:
        i += length
    elif i > length:
        i = length
    return i

def normalize_slice(s, e, length):
    if s is None:
        s = 0
    if e is None:
        e = length
    s = _normalize_slice_h(s, length)
    e = _normalize_slice_h(e, length)

    if s > e:
        s = e = 0
    return range(s, e)

def list_type(L):
    # Return the type of the elements of this list.  If the list is
    # empty or there are multiple types, then return None.
    if not L:
        return None
    # Use the iterator protocol, so it will work with generic
    # sequences (e.g. lists, sets, dicts).
    t = None
    for l in L:
        if t is None:
            t = type(l)
        elif type(l) is not t:
            return None
    return t


## # synonyms should be part of the Matrix class.  Having it separate
## # makes it too hard to create Matrices.
## class _SynonymDecorator:
##     # Members:
##     # _matrix
##     # _synonyms
##     def __init__(self, matrix, synonyms):
##         row_names = {}.fromkeys(matrix.row_names())
##         col_names = {}.fromkeys(matrix.col_names())

##         for synonym, name in synonyms.iteritems():
##             assert name in row_names or name in col_names, "unknown: %s" % name
##         self.__dict__["_matrix"] = matrix
##         self.__dict__["_synonyms"] = synonyms.copy()
##         self.__dict__["_method"] = None
##         #self._matrix = matrix
##         #self._synonyms = synonyms.copy()
##         #self._method = None
##     def _decorated_msv(
##         self, row=None, col=None, row_header=None, col_header=None):
##         if(row_header and row_header not in self._matrix.row_names() and 
##            row_header in self._synonyms):
##             row_header = self._synonyms[row_header]
##         if(col_header and col_header not in self._matrix.col_names() and 
##            col_header in self._synonyms):
##             col_header = self._synonyms[col_header]
##         fn = getattr(self._matrix, self._method)
##         assert fn, "Missing: %s" % self._method
##         x = fn(row=row, col=col, row_header=row_header, col_header=col_header)

##         # If the method creates a new Matrix, then make sure that is
##         # also decorated with the same synonyms.
##         if self._method == "matrix":
##             x = add_synonyms(x, self._synonyms)
##         return x
##     def _decorated_rc_names(self, header=None):
##         if header:
##             names = self._matrix.row_names()
##             if self._method == "col_names":
##                 names = self._matrix.col_names()
##             if header not in names and header in self._synonyms:
##                 header = self._synonyms[header]
##         fn = getattr(self._matrix, self._method)
##         assert fn, "Missing: %s" % self._method
##         x = fn(header=header)
##         return x
##     def __setattr__(self, name, value):
##         setattr(self._matrix, name, value)
##     def __getattr__(self, attr):
##         if attr in ["matrix", "slice", "value", "_index"]:
##             self.__dict__["_method"] = attr
##             #self._method = attr
##             return self._decorated_msv
##         elif attr in ["row_names", "col_names"]:
##             self.__dict__["_method"] = attr
##             #self._method = attr
##             return self._decorated_rc_names
##         return getattr(self._matrix, attr)

## def add_synonyms(matrix, synonyms):
##     # synonyms is a dictionary of synonym -> official name.
    
##     # If matrix is already a _SynonymDecorator, then don't decorate it
##     # twice.
##     if matrix.__class__.__name__ == '_SynonymDecorator':
##         x = matrix._synonyms.copy()
##         x.update(synonyms)
##         synonyms = x
##         matrix = matrix._matrix
##     return _SynonymDecorator(matrix, synonyms)

def is_Matrix(x):
    # Check if it's a Matrix object.  Hard to check based on class
    # name because of derived classes and decorators.  So do a quick
    # and dirty check for some expected "matrix" methods.
    import types
    
    METHODS = [
        "matrix", "slice", "value", "dim", "row_names", "col_names",
        ]
    if type(x) is not types.InstanceType:
        return False
    for m in METHODS:
        if not getattr(x, m, None):
            return False
    return True

def _index_strings(words):
    # Return dictionary of string -> list of indexes.
    word2indexes = {}
    for i, word in enumerate(words):
    #for i in range(len(names)):  # slower
        word = words[i]
        if word not in word2indexes:
            word2indexes[word] = [i]
        else:
            word2indexes[word] += [i]
    return word2indexes

def test_Matrix():
    X = [[0, 1, 2], [3, 4, 5]]
    row_names = {
        "A1" : ["G1", "G2"],
        "A2" : ["G1", "G2"],
        "A3" : ["G", "G"],
        }
    col_names = {
        "A1" : ["S1", "S2", "S3"],
        "J1" : ["J1", "J2", "J3"],
        "F1" : ["F1", "F2", "F3"],
        }
    synonyms = {
        "A1s" : "A1",
        "A1s2" : "A1",
        "J1s" : "J1",
        }
    X = InMemoryMatrix(
        X, row_names=row_names, col_names=col_names, synonyms=synonyms)
        #X, row_names=row_names, col_names=col_names, synonyms={})
    Xs = X
    #Xs = add_synonyms(X, synonyms)

    def d(**keywds):
        return keywds
    
    tests = [
        (X.dim, (), {}, (2, 3)),                 # dim method
        (X.slice, (1, 2), {}, [[5]]),            # simple indexing
        (X.slice, (-1, -2), {}, [[4]]),          # negative indexes
        (X.slice, (1, 3), {}, "exception"),      # out of range
        (X.slice, (-5, 0), {}, "exception"),     # negative out of range
        (X.slice, (0, (0, -1)), {}, [[0, 1]]),   # tuples
        (X.slice, (0, (-1, 0)), {}, []),         # bad tuple
        (X.slice, ((-1, 0), (0, 1)), {}, []),    # tuple with tuple
        (X.slice, ((-5, 5), 0), {}, [[0], [3]]), # tuple out of range
        (X.slice, (None, 0), {}, [[0], [3]]),    # None
        (X.slice, ((None, 1), None), {}, [[0, 1, 2]]),          # None, tuple
        (X.slice, (1, [0, 2]), {}, [[3, 5]]),                   # list
        (X.slice, ([1, 0], [2, 0]), {}, [[5, 3], [2, 0]]),      # list order
        (X.slice, (), d(row="G1",row_header="A2"), [[0, 1, 2]]), # name
        (X.slice, (), d(row="G2"), [[3, 4, 5]]),             # row_header=None
        (X.slice, (), d(row="G2",row_header="A2",col="J2",col_header="J1"),
         [[4]]),                                             # row+col names
        (X.slice, (0,), d(col=["S1","J1","J3"]), [[0, 2]]),  # list of cols
        (X.slice, (), d(row="G1", row_header=["A3"]), []),   # annot non-match
        (Xs.slice, (), d(row="G1", row_header=["A3"]), []),  # annot non-match
        (X.value, (), d(row="G1", row_header="A1", col="S1", col_header="A1"),
         [[0]]),                                        # same row/col names
        (X.slice, (), d(col=["F2", "S1"], col_header=["A1", "F1"]),
         [[1, 0], [4, 3]]),                             # reorder col
        (X.slice, (), d(row="X1",row_header="A2"), []),  # missing id
        (X.slice, (), d(row="G",row_header="A3"),
         [[0, 1, 2], [3, 4, 5]]),                       # dup ids
        (X.slice, (1, 2), d(row="G1",row_header="A2"), "exception"),# bad args
        (X.slice, (0,), d(col="J3",col_header="J1"), [[2]]),        # mixed arg
        (X.value, (1, -2), {}, 4),                    # type; 2 ints
        (X.value, ((0, 5), 1), {}, [1, 4]),           # type; tuple
        (X.__getitem__, ("A1",), {}, "exception"),    # [<str>], row annot
        (X.__getitem__, (("G1", "S1"),), {}, [[0]]),     # [<str>, <str>]
        (X.__getitem__, ((["G1", "G2"], "S1"),), {}, [[0], [3]]),
        (X.__getitem__, (("G", "J2"),), {}, [[1], [4]]),
        #(X.__getitem__, (("ABC", "C2"),), {}, 2),     # [<str>, <str>]
        #(X.__getitem__, ((["ABC", "DEF"], "C0"),), {},
        # [0, 3]),                                     # [<str>, <list>]
        (X.matrix, (0, 0), {}, [[0]]),                # .matrix
        (InMemoryMatrix, ([],), {}, "[]"),            # empty matrix
        (X.matrix, ("Z11", ""), {}, "[]"),            # slice empty matrix

        (Xs.slice, (), d(row="G1",row_header="A1"), [[0, 1, 2]]), # synonym
        (Xs.slice, (), d(row="G1",row_header="A1s"), [[0, 1, 2]]),
        (Xs.slice, (), d(row="G1",row_header="J1s"), "exception"),
        (Xs.slice, (), d(col="J2",col_header="J1s"), [[1], [4]]),
        ]

    for fn, args, keywds, gold_standard in tests:
        try:
            output = fn(*args, **keywds)
        except Exception, x:
            output = "exception"
            if output != gold_standard:
                raise
            #raise
        status = "PASSED"
        if str(output) != str(gold_standard):
            status = "FAILED"
        x = status, args, output, gold_standard
        print "\t".join(map(str, x))

def test__index_strings():
    names = ["a", "b", "c", "c", "hello"]
    print _index_strings(names)
    for i in range(10000):
        print i
        for j in range(1000000):
            _index_strings(names)

try:
    import cMatrix
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in cMatrix.__dict__.keys():
        if name.startswith("__"):
            continue
        this_module.__dict__[name] = cMatrix.__dict__[name]

if __name__ == '__main__':
    test_Matrix()
    #test__index_strings()
