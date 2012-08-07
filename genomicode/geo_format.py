"""

Locator strings are tab-delimited strings with the columns:
- "geo_format"
- geo_path
- GSEID
- GPLID
- filename  (optional)

Classes:
GeoMatrix

Functions:
make_locator
is_format
is_matrix
read
write

"""
import os, sys

import arrayio
from genomicode import Matrix

SAMPLE_NAME = arrayio.tdf.SAMPLE_NAME


class GeoMatrix(Matrix.AbstractMatrix):
    def __init__(self, geo_path, GSEID, GPLID, filename=None,
                 datatype=None, subrows=None, subcols=None,
                 synonyms=None):
        # filename is an optional path to the signal_rank_bin file.
        # subrows and subcols are optional lists of indexes to specify
        # a subset of this matrix.
        from genomicode import filelib
        import arrayio
        import geolib
        
        Matrix.AbstractMatrix.__init__(self, synonyms)
        self.geo_path = geo_path
        self.GSEID, self.GPLID = GSEID, GPLID

        if filename is None:
            filename = geolib.FileFactory.find_signal_rank_file(
                geo_path, GSEID, GPLID)
        assert filename is not None, "Missing signal rank file"
        self.filename = filename
        self.handle = open(filename, "rb")
        ## Test whether we can load the files faster from Redis.
        ## import StringIO
        ## import redis
        ## key = "%s-%s" % (self.GSEID, self.GPLID)
        ## r = redis.Redis(host='localhost', port=6379, db=0)
        ## x = r.get(key)
        ## self.handle = StringIO.StringIO(x)
        
        self.datatype = datatype
        nrow, ncol = geolib.load_signal_bin_matrix_size(self.handle)
        if subrows is None:
            subrows = range(nrow)
        if subcols is None:
            subcols = range(ncol)
        assert not subrows or max(subrows) < nrow
        assert not subcols or max(subcols) < ncol
        self.subrows = list(subrows)
        self.subcols = list(subcols)
        self._num_rows = len(self.subrows)    # optimization
        self._num_cols = len(self.subcols)

        # For optimization.
        self._full_rows = (subrows == range(nrow))
        self._full_cols = (subcols == range(ncol))

        self._rownames = None
        self._rowcache = {}
        self._cache = {}   # index -> expression values for row
        
        # Set the cache to hold at most X numbers in memory at a time.
        # Test on GSE2109-GPL570, 1785 samples.
        #   X  MEMORY
        #  1m   250m
        #  5m   600m
        # 10m   800m
        MAX_NUMBERS = 2E6
        # MAX_CACHE_SIZE is the number of rows.
        self.MAX_CACHE_SIZE = int(max(10, MAX_NUMBERS/ncol))

    def row_names(self, header=None):
        import geolib
        
        if header is None:
            if self._rownames is None:
                self._rownames = geolib.load_signal_bin_rownames(self.handle)
            return self._rownames

        if header not in self._rowcache:
            ids = geolib.load_signal_bin_rownames(self.handle, header)
            if not self._full_rows:
                ids = [ids[i] for i in self.subrows]
            self._rowcache[header] = ids
        return self._rowcache[header]

    def col_names(self, header=None):
        import geolib
        if header is None:
            return [SAMPLE_NAME]
        assert header == SAMPLE_NAME
        ids = geolib.load_signal_bin_sample_ids(self.handle)
        if not self._full_cols:
            ids = [ids[i] for i in self.subcols]
        return ids
    
    def _matrix(self, row_indexes, col_indexes):
        # Return another instance of this object.
        assert not row_indexes or max(row_indexes) < self._num_rows
        assert not col_indexes or max(col_indexes) < self._num_cols
        subrows, subcols = self.subrows, self.subcols
        if row_indexes is not None:
            subrows = [subrows[i] for i in row_indexes]
        if col_indexes is not None:
            subcols = [subcols[i] for i in col_indexes]
        x = GeoMatrix(
            self.geo_path, self.GSEID, self.GPLID, filename=self.filename,
            datatype=self.datatype, subrows=subrows, subcols=subcols)
        x._cache = self._cache  # provide a pointer to share caches.
        x._rownames = self._rownames
        for header, ids in self._rowcache.iteritems():
            if row_indexes is None:
                x._rowcache[header] = ids
            else:
                x._rowcache[header] = [ids[i] for i in row_indexes]
        return x
    
    def _load_rows_cached(self, indexes):
        import geolib
        I = [i for i in indexes if i not in self._cache]
        if I:
            X = geolib.load_signal_bin_many_genes(self.handle, I)
            for i, x in zip(I, X):
                self._cache[i] = x
        X = [self._cache[i] for i in indexes]

        # Clear the cache if it's too big.
        indexes_dict = {}.fromkeys(indexes)
        if len(self._cache) > self.MAX_CACHE_SIZE:
            keys_to_clear = self._cache.keys()[:self.MAX_CACHE_SIZE/2]
            for key in keys_to_clear:
                del self._cache[key]

        return X
    
    def _slice(self, row_indexes, col_indexes):
        # Return just the underlying matrix.
        subrows, subcols = self.subrows, self.subcols
        if row_indexes is not None:
            subrows = [subrows[i] for i in row_indexes]
        if col_indexes is not None:
            subcols = [subcols[i] for i in col_indexes]
        ## Not worth optimizing by taking the subset only when
        ## necessary.
        #subrows = [self.subrows[i] for i in row_indexes]
        #subcols = [self.subcols[i] for i in col_indexes]
        #X = geolib.load_signal_bin_many_genes(self.handle, subrows)
        X = self._load_rows_cached(subrows)
        if X and subcols != range(len(X[0])):
            X = [[x[i] for i in subcols] for x in X]
        return X
    
    def dim(self):
        # Return tuple of (nrow, ncol).
        return self._num_rows, self._num_cols

def make_locator(geo_path, GSEID, GPLID, filename=None):
    x = ["geo_format", geo_path, GSEID, GPLID]
    if filename is not None:
        x.append(filename)
    x = ":".join(x)
    return x

def is_format(locator_str, hrows=None, hcols=None):
    if locator_str.startswith("geo_format"):
        return True
    return False

def is_matrix(X):
    if not hasattr(X, "filename"):
        return False
    if not hasattr(X, "subrows") or not hasattr(X, "subcols"):
        return False
    if not hasattr(X, "MAX_CACHE_SIZE"):
        return False
    return True

def read(handle, hrows=None, hcols=None, datatype=None):
    import arrayio
    
    assert type(handle) is type("")
    assert is_format(handle)
    parts = handle.split(":")
    assert len(parts) in [4, 5]
    x, geo_path, GSEID, GPLID = parts[:4]
    filename = None
    if len(parts) == 5:
        filename = parts[4]
    #print handle, filename

    synonyms = {
        arrayio.ROW_ID : "Probe.Set.ID",
        arrayio.GENE_ID : "LocusLink",
        arrayio.GENE_SYMBOL : "Gene.Symbol",
        arrayio.COL_ID : SAMPLE_NAME,
        }

    X = GeoMatrix(
        geo_path, GSEID, GPLID, filename=filename, datatype=datatype,
        synonyms=synonyms)
    return X

def write(X, handle):
    raise NotImplementedError, "Saving of GEO DataSets not supported"

def _geo_to_tdf(X):
    import arrayio
    from genomicode import Matrix
    
    assert is_matrix(X)
    assert X.row_names()
    assert X.col_names()

    _X = X.slice()
    row_order = X.row_names()
    col_order = X.col_names()
    row_names = {}
    col_names = {}
    synonyms = {}

    for name in X.row_names():
        row_names[name] = X.row_names(name)
    for name in X.col_names():
        col_names[name] = X.col_names(name)

    # Find a suitable ROW_ID.
    row_ids = ["Probe.Set.ID", "Probe Set ID", "PSID", "ID"]
    for row_id in row_ids:
        if row_id in row_order:
            break
    else:
        row_id = row_order[0]
    # Make row_id the first column.
    i = row_order.index(row_id)
    row_order.pop(i)
    row_order = [row_id] + row_order

    # Set the synonyms.
    synonyms[arrayio.ROW_ID] = row_id
    synonyms[arrayio.COL_ID] = col_order[0]
    if "Gene.Symbol" in row_order:
        synonyms[arrayio.GENE_SYMBOL] = "Gene.Symbol"
    if "Description" in row_order:
        synonyms[arrayio.GENE_DESCRIPTION] = "Description"
    if "LocusLink" in row_order:
        synonyms[arrayio.GENE_ID] = "LocusLink"
    
    x = Matrix.InMemoryMatrix(
        _X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #x = Matrix.add_synonyms(x, synonyms)
    assert arrayio.tab_delimited_format.is_matrix(x)
    return x


this_module = sys.modules[__name__]
if this_module not in arrayio.FORMATS:
    arrayio.FORMAT_NAMES.insert(0, sys.modules[__name__])
    arrayio.FORMATS.insert(0, this_module)
    x = "geo_format", "tab_delimited_format", _geo_to_tdf
    arrayio.CONVERTERS.insert(0, x)
del this_module
