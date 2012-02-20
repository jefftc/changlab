"""

Functions:
lwrite       Write to a handle, locking it to prevent concurrent writing.
tswrite      Write to a handle with a timestamp.
openfh       Open a file name or handle.
exists       Whether a filename exists.  Also checks for .gz and .bz2.
exists_nz    Whether a filename exists and has non-zero size.
safe_unlink  Unlink file only if it exists.

read_row     Read one row from a tab-delimited table from a file.
write_row    Save one row from a tab-delimited table to a file.
read_cols

split_by     Split a column.
join_by      Join a column.

iter_by      Group sequential records based on some key value.
as_dict      Convert records into a dictionary.

"""
import os, sys

# _lock
# _unlock
#
# _parse_format
# _read_fmt2fn
# _write_fmt2fn
# _make_convert_fns

class GenericObject:
    def __init__(self, **keywds):
        #import traceback; traceback.print_stack()
        for name, value in keywds.iteritems():
            setattr(self, name, value)
    def __repr__(self):
        # Bug: Does not properly quote strings.
        x = ["  %s=%s" % (n, v) for (n, v) in self.__dict__.iteritems()]
        x = ",\n".join(x)
        return "GenericObject(\n%s\n  )" % x

def _lock(handle):
    import fcntl
    import time

    return # XXX
    fileno = handle.fileno()
    start = time.time()
    while 1:
        try:
            fcntl.lockf(fileno, fcntl.LOCK_EX)
        except Exception, x:
            if str(x).find('No record locks available') < 0 or \
               time.time() >= start + 600:   # try to lock for 10 minutes
                raise
        else:
            break

def _unlock(handle):
    import fcntl

    return # XXX
    fileno = handle.fileno()
    fcntl.lockf(fileno, fcntl.LOCK_UN)

def lwrite(s, handle=None):
    if handle is None:
        import sys
        handle = sys.stdout
    _lock(handle)
    try:
        handle.write(s)
        handle.flush()
    finally:
        _unlock(handle)

def nlwrite(s, handle=None):
    if handle is None:
        import sys
        handle = sys.stdout
    handle.write(s)
    
def tswrite(s, handle=None, format="%m/%d/%Y %H:%M:%S", write_fn=lwrite):
    import time
    time_tup = time.localtime(time.time())
    now = time.strftime(format, time_tup)
    write_fn("%s\t%s" % (now, s), handle=handle)

def _my_popen(cmd):
    import sys
    if sys.hexversion >= 0x02040000:
        from subprocess import Popen, PIPE
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE,
            close_fds=True, universal_newlines=True)
        w, r, e = p.stdin, p.stdout, p.stderr
    else:
        import os
        w, r, e = os.popen3(cmd)
    return w, r, e
        

def openfh(file_or_handle, mode='rU'):
    # Bug: This doesn't handle newlines correctly.  Most of the
    # utilities here will split only on newlines (\n).  So Mac
    # formatted files that have only carriage returns (\r) will not
    # work correctly.  Now fixed.  Uses subprocess module.
    if type(file_or_handle) is not type(''):
        # If this is not a string, assume it's already a file handle.
        return file_or_handle
    elif file_or_handle.lower().startswith("http"):
        # Looks like a URL.
        import urllib2
        return urllib2.urlopen(file_or_handle)
    elif file_or_handle.lower().endswith(".gz"):
        if "r" in mode:
            # May cause broken pipes.  Thus, capture stderr and get
            # rid of it.
            #return os.popen("zcat %s" % file_or_handle)
            if not os.path.exists(file_or_handle):
                raise IOError, "File does not exist: %r" % file_or_handle
            # If file isn't finished reading, should close it so
            # process doesn't stick around half done.
            cmd = "gunzip -c '%s'" % file_or_handle
            w, r, e = _my_popen(cmd)
            w.close()
            e.close()
            return r
        else:
            import gzip
            return gzip.open(file_or_handle, mode)
    elif file_or_handle.lower().endswith(".bz2"):
        if "r" in mode:
            if not os.path.exists(file_or_handle):
                raise IOError, "File does not exist: %s" % file_or_handle
            cmd = "bzcat '%s'" % file_or_handle
            w, r, e = _my_popen(cmd)
            w.close()
            e.close()
            return r
        else:
            raise NotImplementedError
    elif file_or_handle.lower().endswith(".zip"):
        if "r" in mode:
            if not os.path.exists(file_or_handle):
                raise IOError, "File does not exist: %s" % file_or_handle
            cmd = "unzip -p '%s'" % file_or_handle
            w, r, e = _my_popen(cmd)
            w.close()
            e.close()
            return r
        else:
            raise NotImplementedError
    elif file_or_handle.lower().endswith(".xls") or \
        file_or_handle.lower().endswith(".xlsx"):
        cmd = "xls2txt '%s'" % file_or_handle
        w, r, e = _my_popen(cmd)
        # This may block.  Just ignore this until we can figure out
        # how to fix it properly.
        #x = e.read()
        #assert not x, x
        w.close()
        e.close()
        return r
    return open(file_or_handle, mode)

def exists(filename):
    if type(filename) is type("") and filename.lower().startswith("http"):
        # Looks like a URL.
        import urllib2
        try:
            handle = urllib2.urlopen(filename)
        except urllib2.HTTPError:
            return None
        return filename

    EXTS = [".gz", ".bz2", ".zip"]
    name, ext = os.path.splitext(filename)
    if ext.lower() in EXTS:
        filename = name
    if os.path.exists(filename):
        return filename
    for ext in EXTS:
        f = filename + ext
        if os.path.exists(f):
            return f
    return None

def exists_nz(filename):
    import stat
    fn = exists(filename)
    if not fn:
        return None
    if os.stat(fn)[stat.ST_SIZE] > 0:
        return fn
    return None
    
def _parse_format(format):
    """Return names, format."""
    assert format is not None
    if format.find(":") < 0:
        assert " " not in format
        return None, format
    
    names, fmt = [], ""
    for x in format.split():
        name, type = x.split(":")
        assert len(type) == 1, "Invalid type: %s" % type
        names.append(name)
        fmt += type
    assert len(names) == len(fmt)
    return names, fmt

def _read_fmt2fn(f):
    # Return a function that can convert this format, or None.
    if f == 's':
        return None
    elif f in ['d', 'i']:
        return int
    elif f == 'l':
        return long
    elif f == 'f':
        return float
    elif f == 'x':
        return None
    raise ValueError, "Unknown format specifier %s" % f

def _write_fmt2fn(f):
    # Return a function that can convert this format, or None.
    if f == 's':
        return str
    elif f in ['d', 'i']:
        return str
    elif f == 'l':
        return str
    elif f == 'f':
        return str
    elif f == 'x':
        return None
    raise ValueError, "Unknown format specifier %s" % f

def _make_convert_fns(format, obj_convert_fns, fmt2fn_fn):
    obj_convert_fns = list(obj_convert_fns)
    convert_fns = []
    for f in format:
        if f == "O":
            fn = obj_convert_fns.pop(0)
        else:
            fn = fmt2fn_fn(f)
        convert_fns.append(fn)
    return convert_fns

# Bug: if the file is gzip'd, will leave gunzip -c processes lying
# around.
def read_cols(file_or_handle, delimiter="\t", skip=0):
    handle = openfh(file_or_handle)
    for i in range(skip):
        handle.readline()
    for line in handle:
        cols = line.rstrip("\r\n").split(delimiter)
        yield cols
    handle.close()

def _normalize_colname(name):
    import re
    # Fix the header to be a python variable.
    x = str(name)
    # Replace all non-word character with _.
    x = re.sub(r"\W", "_", x)
    # Replace initial numbers with Xnumber.
    x = re.sub(r"^(\d)", r"X\1", x)
    return x

def _make_format_from_header(names):
    import math
    
    # Normalize each name.
    normnames = [_normalize_colname(x) for x in names]

    # For columns with duplicate names, number them so that they are
    # unique.
    name2count = {}
    for n in normnames:
        name2count[n] = name2count.get(n, 0) + 1
    name2nextid = {}
    for i, name in enumerate(normnames):
        count = name2count[name]
        if count == 1:
            continue
        # Figure out the number of digits to use for the id.
        ndigits = int(math.floor(math.log(count, 10))) + 1
        id = name2nextid.get(name, 0)
        name2nextid[name] = id + 1
        x = "%s_%0*d" % (name, ndigits, id)
        normnames[i] = x

    # Make each of the columns a string.
    normnames = ["%s:s" % x for x in normnames]

    # Return the format.
    return " ".join(normnames)
                
class RowIterator:
    # Member variables:
    # _line        Previous line from the file (unparsed).
    # _cols
    # _header      Name of each column.  From header, format, then index.
    # _nheader     Normalized names.

    def __init__(self, file_or_handle, delimiter, strip, skip, header, 
                 format=None, *obj_convert_fns):
        if skip:
            file_or_handle = openfh(file_or_handle)
            for i in range(skip):
                file_or_handle.readline()
        reader = self._parse_line(file_or_handle, delimiter=delimiter)

        names = None
        if header:
            # If the file is empty, this will raise a StopIteratio
            # exception.
            try:
                line, names = reader.next()
            except StopIteration:
                names = []
            # If no format is provided, then make one from the header.
            if format is None:
                format = _make_format_from_header(names)
        assert format is not None, "No format given."

        normnames, format = _parse_format(format)
        if not names:
            names = normnames
        x = _make_convert_fns(format, obj_convert_fns, _read_fmt2fn)
        convert_fns = x
        fn_cols = [
            i for (i, fn, f)
            in zip(range(len(format)), convert_fns, format)
            if fn != None and f != "x"]
        x = [i for (i, f) in enumerate(format) if strip and f == "s"]
        strip_cols = x

        self._header = names
        self._nheader = normnames
        self._reader = reader
        self._format = format
        self._convert_fns = convert_fns
        self._fn_cols = fn_cols
        self._strip_cols = strip_cols

    def _parse_line(self, file, delimiter="\t"):
        import csv
        # Allow up to 32Mb fields (Python 2.5 and above).
        if hasattr(csv, "field_size_limit"):
            csv.field_size_limit(32*1024*1024)
        for i, line in enumerate(openfh(file)):
            # Skip blank lines.
            if not line.strip():
                continue
            try:
                x = csv.reader([line], delimiter=delimiter).next()
            except csv.Error, x:
                raise csv.Error, "%s [%d: %s]" % (str(x), i, repr(line))
            yield line, x

    def next(self):
        line, data = self._reader.next()
        cols = data[:]

        if len(data) != len(self._format):
            dlen, flen = len(data), len(self._format)
            s = "data(%d)/format(%d) are different lengths\n%r\n%r" % (
                dlen, flen, self._format, data)
            raise AssertionError, s
        for i in self._fn_cols:
            data[i] = self._convert_fns[i](data[i])
        for i in self._strip_cols:
            data[i] = data[i].strip()
        if self._nheader:
            params = {}
            for (n, d, f) in zip(self._nheader, data, self._format):
                if f == "x":
                    continue
                params[n] = d
            data = GenericObject(**params)

        self._line = line
        self._cols = cols
        setattr(data, "_iter", self)
        setattr(data, "_line", line)
        setattr(data, "_cols", cols)
        setattr(data, "_header", self._header)
        setattr(data, "_nheader", self._nheader)
        return data

    def __iter__(self):
        return self

# read_row(filename, "ssss")
# read_row(filename, "name:s number:i")
# read_row(filename, "col1:O col2:O", convert_fn, convert_fn)
# read_row(filename, header=1)
#   delimiter   Character used for the delimiter.  default "\t".
#   skip        Number of lines to skip at the beginning.
def read_row(file_or_handle, *args, **keywds):
    """Iterate over each line of a tab-delimited file.  The iterator
    object contains the following member variables:
        _line       Previous line from the file (unparsed).
        _cols       Previous line from the file, parsed as columns.
        _header     List of column names.  From header, format, then index.
        _nheader    Normalized names.
        
    format is a string describing the type of each column.  The format
    string can follow one of two syntaxes:
    
    1.  "<type1><type2><type3> ..."
    2.  "<name1>:<type1> <name2>:<type2> ..."

    In the first syntax, each character in the string corresponds to
    the type of a column.  There should be no intervening spaces.  The
    second syntax also corresponds to columns, but names are also
    given.

    If names are given, then returns an object with the names used as
    member variables.  Otherwise, returns a list.  Either way, the
    returned object will contain special member variables:
      _iter   reference back to the original iterator
      _line
      _cols
      _header
      _nheader
    
    Allowed types are:
    s  string
    i  integer
    f  float
    O  object
    x  ignore this column

    Each "object" type requires a corresponding function to convert.

    If header=1 is given, will use it for the format string.  If the
    format string is provided, then will just skip the header.

    """
    known_params = ["delimiter", "strip", "skip", "header"]
    for key in keywds:
        assert key in known_params, "Unknown parameter: %s" % key
    delimiter = keywds.get("delimiter", "\t")
    strip = keywds.get("strip", False)
    skip = keywds.get("skip", 0)
    header = keywds.get("header", False)
    return RowIterator(file_or_handle, delimiter, strip, skip, header, *args)

def write_row(file_or_handle, data, format=None, *obj_convert_fns):
    """Write one row of a table into file_or_handle.  data should be a
    list or tuple.  It can be a dict, if the format contains names."""
    import csv
    import operator

    if file_or_handle is None:
        file_or_handle = sys.stdout
    assert type(file_or_handle) is not type(""), "Can not write row to a file."
    # If named_format is given, should I accept a dict for data?
    names, format = _parse_format(format)
    if not operator.isSequenceType(data):  # assume it's a TableRow
        assert names
        data = [getattr(data, n) for n in names]
    assert len(data) == len(format), "data/format are different lengths"
    convert_fns = _make_convert_fns(format, obj_convert_fns, _write_fmt2fn)
    
    handle = openfh(file_or_handle, 'w')
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    # Should ignore the columns with "x" format, so get rid of it.
    data = [fn(d) for (d, f, fn) in zip(data, format, convert_fns) if f != 'x']
    w.writerow(data)
    #handle.close()

def split_by(delimiter, convert_fn=None):
    def f(s):
        x = s.split(delimiter)
        if convert_fn is not None:
            x = [convert_fn(x) for x in x]
        return x
    return f

def join_by(delimiter, convert_fn=None):
    def f(data):
        if convert_fn is not None:
            data = [convert_fn(x) for x in data]
        return ";".join(data)
    return f

def iter_by(iterator, *args, **params):
    assert args, "no fieldnames provided"
    VALID_PARAMS = ["BATCH_SIZE"]
    for k in params:
        assert k in VALID_PARAMS
    BATCH_SIZE = params.get("BATCH_SIZE", None)
    
    started = 0
    prev_key, data = None, []
    for d in iterator:
        key = [getattr(d, fn) for fn in args]
        if (not started or key != prev_key or
            (BATCH_SIZE is not None and len(data) >= BATCH_SIZE)):
            if data:
                yield data
            started = 1
            prev_key, data = key, []
        data.append(d)
    if data:
        yield data

def as_dict(iterator, *args, **keywds):
    # args is the member variables that comprise the key.
    # 
    # Keyword arguments:
    # unique     Save the values as an object (not a list).
    # use_first  Only save the first value for each key.
    unique = keywds.get("unique", 1)
    use_first = keywds.get("use_first", 0)

    dict = {}
    for d in iterator:
        key = tuple([getattr(d, fn) for fn in args])
        if len(args) == 1:
            key = key[0]
        if not unique:
            dict.setdefault(key, []).append(d)
        elif not (use_first and key in dict):
            dict[key] = d
    return dict

def safe_unlink(filename):
    if not filename or not os.path.exists(filename):
        return
    os.unlink(filename)

class DelFile:
    def __init__(self, filename, mode):
        self.filename = filename
        self.handle = open(filename, mode)
    def __getattr__(self, attr):
        return getattr(self.handle, attr)
    def __del__(self):
        if not self.handle.closed:
            self.handle.close()
        safe_unlink(self.filename)

def make_temp_handle(suffix="", prefix="", dir=None, realfile=True):
    from StringIO import StringIO
    import tempfile

    if realfile:
        x, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir)
        handle = DelFile(filename, "r+w")
    else:
        handle = StringIO()
    return handle
