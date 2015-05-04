"""

Objects:
FileFactory

Functions:
download_seriesmatrix_file

align
filter_mv

save_signal_bin
load_signal_bin_matrix_size
load_signal_bin_sample_ids
load_signal_bin_rownames
load_signal_bin_gene
load_signal_bin_many_genes

load_alpha_z_cutoff

read_and_unpack
pack_and_write
read_and_unpack_str
pack_and_write_str

"""
# Internal Functions:
# _extract_sm_sample_meta
# _clean_sm_sample_meta
# _prettify_sm_sample_meta
# 
# _remove_quotes

import os, sys

class FileFactory:
    def __init__(self):
        self._path = self._make_fn = None
        self._extension = None
        
    def _make_file(self, path, name, force_gzip_ext):
        assert self._path
        p = os.path.join(path, self._path)
        if not os.path.exists(p):
            os.mkdir(p)
        filename = os.path.join(p, name)
        if force_gzip_ext:
            filename = filename + ".gz"
        self._path = None
        return filename

    def _check_and_normalize_id(self, id, prefix):
        id = id.upper()
        assert id.startswith(prefix), id
        return id

    def _make_gseid_file(self, path, gseid, force_gzip_ext=0):
        assert self._extension
        gseid = self._check_and_normalize_id(gseid, "GSE")
        return self._make_file(
            path, "%s.%s" % (gseid, self._extension), force_gzip_ext)
    def _make_gplid_file(self, path, gplid, force_gzip_ext=0):
        assert self._extension
        gplid = self._check_and_normalize_id(gplid, "GPL")
        return self._make_file(
            path, "%s.%s" % (gplid, self._extension), force_gzip_ext)
    def _make_gseid_gplid_file(self, path, gseid, gplid, force_gzip_ext=0):
        assert self._extension
        gseid = self._check_and_normalize_id(gseid, "GSE")
        gplid = self._check_and_normalize_id(gplid, "GPL")
        return self._make_file(
            path, "%s-%s.%s" % (gseid, gplid, self._extension), force_gzip_ext)

    def _find_file(self, *attrs, **keywds):
        # Return filename or None.
        import stat
        
        assert self._make_fn
        filename = self._make_fn(*attrs, **keywds)
        self._make_fn = None
        exts = ["", ".gz", ".bz2"]
        for ext in exts:
            if os.path.exists(filename+ext) and \
                   os.stat(filename+ext)[stat.ST_SIZE] > 0:
                break
        else:
            return None
        return filename+ext
    
    def __getattr__(self, attr):
        # make_<name>_file
        # find_<name>_file
        assert attr.startswith("make_") or attr.startswith("find_")
        assert attr.endswith("_file")
        name = attr[5:-len("_file")]
        name2data = {
            "family_soft" : ("family_soft", 1, 0, "txt"),
            "platform" : ("platform", 0, 1, "txt"),
            "annot" : ("annot", 0, 1, "txt"),
            "series_matrix" : ("series_matrix", 1, 1, "txt"),

            "signal" : ("signal_txt", 1, 1, "txt"),
            "signal_bin" : ("signal_bin", 1, 1, "bin"),
            "signal_rank" : ("signal_rank", 1, 1, "bin"),
            
            #"signal" : ("signal_raw", 1, 1),
            #"signal_bin" : ("signal_raw_bin", 1, 1),
            #"signal_rank" : ("signal_rank", 1, 1),
            #"signal_rank_bin" : ("signal_rank_bin", 1, 1),
            #"signal_annot" : ("signal_annot", 1, 1),
            #"signal_annot_bin" : ("signal_annot_bin", 1, 1),
            #"annot_indexes" : ("annot_indexes", 1, 1),
            
            "shapiro" : ("shapiro", 1, 1, "txt"),
            "gene_rank" : ("gene_rank", 1, 1, "txt"),
            "alpha_dist" : ("alpha_dist", 1, 1, "txt"),
            "signal_hg18" : ("signal_hg18", 1, 1, "bin"),
            }
        path, use_gseid, use_gplid, extension = name2data[name]
        self._extension = extension
        if use_gseid and use_gplid:
            fn = self._make_gseid_gplid_file
        elif use_gseid:
            fn = self._make_gseid_file
        elif use_gplid:
            fn = self._make_gplid_file
        else:
            raise AssertionError
        self._path = path
        if attr.startswith("find_"):
            self._make_fn = fn
            fn = self._find_file
        return fn

FileFactory = FileFactory()

## def load_filtered_signals(
##     geo_path, gseid, gplid, percent_filter_by_mean, percent_filter_by_var):
##     import matrixfns as mf

##     filename = FileFactory.find_signal_file(geo_path, gseid, gplid)
##     signal_data = mf.load_arrayfile(
##         filename, num_header_rows=1, num_header_cols=1)
    
##     filename = FileFactory.find_gene_rank_file(geo_path, gseid, gplid)
##     rank_data = mf.load_arrayfile(
##         filename, num_header_rows=1, num_header_cols=1)

##     signal_data = filter_mv(
##         signal_data, rank_data, percent_filter_by_mean, percent_filter_by_var)
##     return signal_data


def download_seriesmatrix_file(outhandle, GSEID, GPLID=None):
    import tempfile
    import gzip
    from ftplib import FTP
    
    assert GSEID.startswith('GSE'), "Bad GSEID: %s" % GSEID
    assert not GPLID or GPLID.startswith('GPL'), "Bad GPLID: %s" % GPLID

    # SeriesMatrix files stored in:
    # ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/<GSEID>/
    #   <GSEID>_series_matrix.txt.gz
    #   <GSEID>-<GPLID>_series_matrix.txt.gz  (if multiple GPLID)

    # Example GSEID
    # GSE9691  Only 1 GPLID
    # GSE6532  Multiple GPLID
    
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()
    ftp.cwd('pub/geo/DATA/SeriesMatrix/' + GSEID)

    # Look for the file in the series matrix path.
    files = []
    ftp.retrlines('NLST', files.append)

    desired_filename = "%s_series_matrix.txt.gz" % GSEID
    if GPLID:
        desired_filename = "%s-%s_series_matrix.txt.gz" % (GSEID, GPLID)
    assert desired_filename in files, (
        "I could not find the series matrix file for %s.  "
        "Possibilities are:\n%s." % (GSEID, "\n".join(files)))

    # Download the gzip'd file.
    gzipped_file = None
    try:
        x, gzipped_file = tempfile.mkstemp(dir=".", suffix=".gz"); os.close(x)
        handle = open(gzipped_file, 'wb')
        ftp.retrbinary('RETR ' + desired_filename, handle.write)
        handle.close()

        # Un-gzip the file.
        handle = gzip.GzipFile(gzipped_file, 'rb')
        while True:
            line = handle.readline()
            if not line:
                break
            outhandle.write(line)
    finally:
        if gzipped_file and os.path.exists(gzipped_file):
            os.unlink(gzipped_file)


def _remove_quotes(s):
    if s and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s


def _extract_sm_sample_meta(file_or_handle):
    # Return a tab-delimited matrix of the results.
    import gzip
    import shlex
    import filelib

    # Sometimes there are \r embedded within quotes.  Split lines
    # based on \n.
    if type(file_or_handle) is type(""):
        filename = file_or_handle
        if filename.endswith(".gz"):
            handle = gzip.open(filename, 'rb')
        elif filename.endswith(".bz2"):
            cmd = "bzcat '%s'" % filename
            w, r, e = filelib._my_popen(cmd)
            w.close()
            e.close()
            handle = r
        else:
            # Assume text file.
            handle = open(filename)
    else:
        handle = file_or_handle
    lines = handle.read().split("\n")

    # Sometimes if a line is too long, it will be continued on the
    # next line, e.g.
    # !Sample_extract_protocol_ch1  <TEXT>  <TEXT>  <TEXT>
    # !Sample_extract_protocol_ch1  <CONT>  <CONT>  <CONT>
    # !Sample_extract_protocol_ch1          <CONT>  <CONT>
    # !Sample_extract_protocol_ch1          <CONT>  <CONT>
    #
    # Notes:
    # o It seems like the !Sample_<name> is consistent.
    # o Some of the columns can be blank if there is no continuation.
    # o Sometimes the close quotes are missing from <TEXT> or <CONT>.
    matrix = []
    num_cols = None
    for line in lines:
        line = line.rstrip("\n")

        if not line.startswith("!Sample"):
            continue

        # Split by \t.
        # Use shlex because sometimes there are tabs within the
        # quotes.  However, this will fail if the close quote is
        # missing.  If this is the case, parse without shlex and hope
        # there aren't any internal tabs.
        if not line.endswith('"'):
            cols = line.split("\t")
        else:
            x = shlex.shlex(line)
            x.whitespace = "\t"
            x.whitespace_split = True
            cols = [x for x in x]
        # Sometimes shlex will mess up if there are quotes within
        # strings.  If this happens, fall back on just splitting on
        # tabs.
        if num_cols and len(cols) != num_cols:
            cols = line.split("\t")
        if line.startswith("!Sample_title"):
            num_cols = None
        if num_cols is None:
            num_cols = len(cols)
        assert len(cols) == num_cols, "%s: has %d expect %d" % (
            cols[0], len(cols), num_cols)
        cols = [_remove_quotes(x) for x in cols]
        cols = [x.replace("\t", " ").strip() for x in cols]
        cols = [x.replace("\r", " ").strip() for x in cols]
        matrix.append(cols)
    return matrix


SAMPLE_CHARACTERISTICS = "!Sample_characteristics_ch1"
def _clean_sm_sample_meta(matrix):
    # Return a cleaned up matrix.  Does two things:
    # 1.  Merge into one table where each column is a single sample.
    #     Sometimes wide matrices (too many columns) get split up.
    #     Merge them back together.
    # 2.  Sample_characteristics_ch1 are not guaranteed to be aligned.
    #     Align them.
    

    # Examples:
    # GSE25066   Alignment messed up.

    #IGNORE = [
    #    "Sample_supplementary_file",
    #    "Sample_description",
    #    "Sample_extract_protocol_ch1",
    #    "Sample_growth_protocol_ch1",
    #    "Sample_status",
    #    "Sample_submission_date",
    #    "Sample_last_update_date",
    #    "Sample_treatment_protocol_ch1",
    #    "Sample_data_row_count",
    #    "Sample_scan_protocol",
    #    ]

    # Remove quotes around elements.
    #for i in range(len(matrix)):
    #    for j in range(len(matrix[i])):
    #        x = matrix[i][j]
    #        if x.startswith('"') and x.endswith('"'):
    #            x = x[1:-1]
    #            matrix[i][j] = x

    # In this sample matrix, each row can only have 256 columns.  If
    # it is longer, it will be broken up across the file.  Fix this.
    while True:
        assert matrix
        assert matrix[0][0] == "!Sample_title", matrix[0][0]
        i_next = None
        for i in range(1, len(matrix)):
            #print i, matrix[i][0]
            if matrix[i][0] == "!Sample_title":
                i_next = i
                break
        if not i_next:  # no more matrices
            break
        
        # Merge this with the previous matrix.
        num_rows = i_next
        for i in range(num_rows):
            # Make sure the header is the same.
            j = i+num_rows
            #print len(matrix), i_next, num_rows, i, j
            assert matrix[i][0] == matrix[j][0]
            x = matrix[i] + matrix[j][1:]
            matrix[i] = x
        # Delete the rows I just merged.
        for i in range(num_rows):
            del matrix[i_next]

    # The !Sample_characteristics_ch1 lines sometimes get unaligned:
    # !Sample_characteristics_ch1  subtype: type VI  time_to_metastasis: 1.6
    # !Sample_characteristics_ch1                    subtype: type IV
    # If this happens, fix the alignment.

    # First, make a list of all the Sample_characteristics_ch1 lines.
    I_char = []
    for i in range(len(matrix)):
        if matrix[i][0] == SAMPLE_CHARACTERISTICS:
            I_char.append(i)
            
    # For each sample, parse out all the sample characteristics.
    # Also, make an in-order list of all the possible
    # sample_characteristics.
    sample2chars = {}  # index -> dict of characteristics
    all_chars = []
    for j in range(1, len(matrix[0])):
        chars = {}
        for i in I_char:
            # Sample_characteristics_ch1 have the format:
            # <name>: <value>
            x = matrix[i][j]
            if not x:
                continue
            assert ":" in x, x
            (name, value) = x.split(":", 1)
            name, value = name.strip(), value.strip()
            chars[name] = value
            if name not in all_chars:
                all_chars.append(name)
        sample2chars[j] = chars

    # Now create the aligned Sample_characteristics lines.
    cmatrix = []
    for i in range(len(all_chars)):
        x = [SAMPLE_CHARACTERISTICS] + [""] * len(sample2chars)
        cmatrix.append(x)
    for j in range(1, len(matrix[0])):
        for i, name in enumerate(all_chars):
            value = sample2chars[j].get(name)
            x = ""
            if value is not None:
                x = "%s: %s" % (name, value)
            cmatrix[i][j] = x
    assert len(cmatrix[0]) == len(matrix[0])

    # Now remove the SAMPLE_CHARACTERISTICS lines from the original
    # matrix.
    omatrix = [x for x in matrix if x[0] != SAMPLE_CHARACTERISTICS]
    assert len(omatrix) >= min(I_char)
    m1 = omatrix[:min(I_char)]
    m2 = omatrix[min(I_char):]
    matrix = m1 + cmatrix + m2

    return matrix


def _prettify_sm_sample_meta(matrix):
    # Return a pretty version of the matrix.
    # 1.  Remove "!Sample_" from the header names.
    #     !Sample_geo_accession -> geo_accession
    # 2.  If only "_ch1" and no "_ch2" suffixes, remove "_ch1".
    # 3.  Rename Sample_characteristics_ch1 to the name given by the
    #     depositor.

    # Do some basic checking on the format.
    assert matrix
    assert matrix[0]
    for i in range(len(matrix)):
        assert matrix[i][0].startswith("!Sample_")

    # Make a copy of the matrix.
    matrix = [x[:] for x in matrix]

    # Remove !Sample_ prefix.
    for i in range(len(matrix)):
        x = matrix[i][0]
        if x == SAMPLE_CHARACTERISTICS:
            continue
        assert x.startswith("!Sample_")
        x = x[len("!Sample_"):]
        assert x
        matrix[i][0] = x

    # Remove _ch1, if there are no _ch2.
    has_ch2 = False
    for i in range(len(matrix)):
        x = matrix[i][0]
        if x.endswith("_ch2"):
            has_ch2 = True
    for i in range(len(matrix)):
        x = matrix[i][0]
        if x == SAMPLE_CHARACTERISTICS:
            continue
        if x.endswith("_ch1") and not has_ch2:
            x = x[:-4]
            assert x
            matrix[i][0] = x

    # Rename SAMPLE_CHARACTERISTICS.
    for i in range(len(matrix)):
        x = matrix[i][0]
        if x != SAMPLE_CHARACTERISTICS:
            continue
        # Figure out the name of this sample characteristic.
        header = None
        for j in range(1, len(matrix[i])):
            x = matrix[i][j]
            if not x:
                continue
            (name, value) = x.split(":", 1)
            name, value = name.strip(), value.strip()
            if header is None:
                header = name
            assert header == name, "unaligned %s line" % SAMPLE_CHARACTERISTICS
            matrix[i][j] = value
        assert header
        matrix[i][0] = header
    
    return matrix
            
        
def align(signal_data, other_data):
    # Align by the IDs in the first column.
    from genomicode import jmath
    assert jmath.nrow(signal_data) <= jmath.nrow(other_data)
    ids1 = jmath.slice(signal_data, (1, jmath.nrow(signal_data)), 0)
    ids2 = jmath.slice(other_data, (1, jmath.nrow(other_data)), 0)
    I = jmath.match(ids1, ids2)
    assert None not in I, "mismatched IDs"
    I = [0] + [x+1 for x in I]   # Take the header too.
    other_data = jmath.slice(other_data, I, None)
    return other_data

def filter_mv(signal_data, rank_data,
              percent_filter_by_mean, percent_filter_by_var):
    from genomicode import jmath
    
    # Align the signal_data to the rank_data.
    rank_data = align(signal_data, rank_data)

    mean_r = jmath.slice(
        rank_data, (1, jmath.nrow(rank_data)), rank_data[0].index("Mean.Rank"))
    var_r = jmath.slice(
        rank_data, (1, jmath.nrow(rank_data)), rank_data[0].index("Var.Rank"))

    mean_n = (1-percent_filter_by_mean) * jmath.nrow(signal_data)
    var_n = (1-percent_filter_by_var) * jmath.nrow(signal_data)
    I = [i for i, (m, v) in enumerate(zip(mean_r, var_r)) if
         m < mean_n and v < var_n]
    raise AssertionError, "Is there an off-by-1 error here?"
    return jmath.slice(signal_data, I, None)

def save_signal_bin(handle, annotations, sample_ids, X):
    # Format:
    # int               num genes
    # int               num samples
    # int               num headers
    # char              typecode of data
    # int               file offset for sample names
    # int               file offset for header names
    # int[num headers]  file offset for annotations i
    # int               file offset for signal data
    #
    # SAMPLE NAMES
    # int               length of sample names
    # char[]            sample names, tab-delimited
    #
    # HEADER NAMES
    # int             length of header names
    # char[]          header names, tab-delimited
    #
    # ANNOTATIONS (REPEAT nheaders TIMES)
    # int             length of annotations
    # char[]          annotations, tab-delimited
    #
    # SIGNAL DATA
    # float[]         gene i sample j x NUMBER OF SAMPLES
    #                 ... x NUMBER OF GENES
    import struct
    from genomicode import jmath
    
    num_genes, num_samples = len(X), len(sample_ids)
    assert num_genes == jmath.nrow(X)
    assert num_samples == jmath.ncol(X)
    for annots in annotations.itervalues():
        assert len(annots) == num_genes
    
    for x in sample_ids:
        assert "\t" not in x
    sample_str = "\t".join(sample_ids)

    # TODO: Should allow user to provide the order of the headers.
    header_names = sorted(annotations)
    assert header_names
    num_headers = len(header_names)
    for x in header_names:
        assert "\t" not in x
    header_str = "\t".join(header_names)

    for name in header_names:
        for x in annotations[name]:
            assert "\t" not in x
    annotation_strs = ["\t".join(annotations[x]) for x in header_names]
    
    sample_offset = (
        struct.calcsize(">2i") +    # num_genes, num_samples
        struct.calcsize(">i") +     # num_headers
        struct.calcsize(">c") +     # typecode
        struct.calcsize(">i") +     # offset for sample names
        struct.calcsize(">i") +     # offset for header names
        struct.calcsize(">%di" % num_headers) +  # offset for annotations
        struct.calcsize(">i")       # offset for signal data
        )
    header_offset = sample_offset + struct.calcsize(">i") + len(sample_str)
    annotation_offsets = [None] * num_headers
    offset = header_offset + struct.calcsize(">i") + len(header_str)
    for i in range(num_headers):
        annotation_offsets[i] = offset
        offset += struct.calcsize(">i") + len(annotation_strs[i])
    signal_offset = offset
    
    # Choose the right typecode.
    typecode = "f"   # float, by default
    if type(X[0][0]) is type(0):
        typecode = "H"      # unsigned short, 2 bytes
        mv = max([max(x) for x in X])
        if mv >= 2**16:
            typecode = "I"  # unsigned int, 4 bytes
        assert mv < 2**32
    
    # Write the header of sizes and offsets.
    pack_and_write(handle, ">2i", num_genes, num_samples)
    pack_and_write(handle, ">i", num_headers)
    pack_and_write(handle, ">c", typecode)
    pack_and_write(handle, ">i", sample_offset)
    pack_and_write(handle, ">i", header_offset)
    pack_and_write(handle, ">%di" % num_headers, *annotation_offsets)
    pack_and_write(handle, ">i", signal_offset)

    pack_and_write_str(handle, sample_str)   # Write the sample names.
    pack_and_write_str(handle, header_str)   # Write the header names.
    for annotation_str in annotation_strs:   # Write the annotations.
        pack_and_write_str(handle, annotation_str)

    # Write the signal data, one gene at a time.
    for x in X:
        x = struct.pack(">%d%s" % (len(x), typecode), *x)
        handle.write(x)

    handle.flush()

def load_signal_bin_matrix_size(handle):
    if type(handle) is type(""):
        handle = open(handle)
    handle.seek(0)
    x = read_and_unpack(handle, ">2i")
    num_genes, num_samples = x
    return num_genes, num_samples

def load_signal_bin_sample_ids(handle):
    if type(handle) is type(""):
        handle = open(handle)

    # Read the header.
    handle.seek(0)
    x = read_and_unpack(handle, ">3ic")
    num_genes, num_samples, num_headers, typecode = x
    x = read_and_unpack(handle, ">ii%dii" % num_headers)
    sample_offset, header_offset = x[:2]
    annotation_offsets = x[2:2+num_headers]
    signal_offset = x[-1]
    #print num_genes, num_samples, num_headers, typecode
    #print sample_offset, header_offset, annotation_offsets, signal_offset
    
    handle.seek(sample_offset)
    x = read_and_unpack_str(handle)
    x = x.split("\t")
    assert len(x) == num_samples
    return x
    
def load_signal_bin_rownames(handle, name=None):
    #print "NAME", name
    if type(handle) is type(""):
        handle = open(handle)

    # Read the header.
    handle.seek(0)
    x = read_and_unpack(handle, ">3ic")
    num_genes, num_samples, num_headers, typecode = x
    x = read_and_unpack(handle, ">ii%dii" % num_headers)
    sample_offset, header_offset = x[:2]
    annotation_offsets = x[2:2+num_headers]
    signal_offset = x[-1]
    assert len(annotation_offsets) == num_headers
    #print num_genes, num_samples, num_headers, typecode
    #print sample_offset, header_offset, annotation_offsets, signal_offset

    # Read the header names.
    handle.seek(header_offset)
    x = read_and_unpack_str(handle)
    header_names = x.split("\t")
    assert len(header_names) == num_headers
    if name is None:
        return header_names
    
    if name not in header_names:
        raise KeyError, "I could not find header: %s" % name
    i = header_names.index(name)
    assert i >= 0 and i < num_headers

    #print "HERE", annotation_offsets
    #print i
    handle.seek(annotation_offsets[i])
    x = read_and_unpack_str(handle)
    x = x.split("\t")
    assert len(x) == num_genes
    return x

def load_signal_bin_gene(handle, index):
    import struct
    if type(handle) is type(""):
        handle = open(handle)
    
    # Read the header.
    handle.seek(0)
    x = read_and_unpack(handle, ">3ic")
    num_genes, num_samples, num_headers, typecode = x
    x = read_and_unpack(handle, ">ii%dii" % num_headers)
    sample_offset, header_offset = x[:2]
    annotation_offsets = x[2:2+num_headers]
    signal_offset = x[-1]

    assert index >= 0 and index < num_genes, "index out of range (%d:%d)" % (
        index, num_genes)

    row_length = struct.calcsize(">%d%s" % (num_samples, typecode))
    offset = signal_offset + row_length*index
    handle.seek(offset)
    return read_and_unpack(handle, ">%d%s" % (num_samples, typecode))

def load_signal_bin_many_genes(handle, indexes):
    import struct
    if type(handle) is type(""):
        handle = open(handle)
        
    # Unfortunately, indexes are scattered all over the file, so have
    # to do a lot of seeks/reads.

    # Read the header.
    handle.seek(0)
    x = read_and_unpack(handle, ">3ic")
    num_genes, num_samples, num_headers, typecode = x
    x = read_and_unpack(handle, ">ii%dii" % num_headers)
    sample_offset, header_offset = x[:2]
    annotation_offsets = x[2:2+num_headers]
    signal_offset = x[-1]

    row_format = ">%d%s" % (num_samples, typecode)
    row_length = struct.calcsize(row_format)

    # Read all the data in to one big string.
    len_indexes = len(indexes)
    raw_data = [None] * len_indexes
    for zzz, i in enumerate(indexes):
        assert i >= 0 and i < num_genes, "index out of range (%d:%d)" % (
            i, num_genes)
        offset = signal_offset + row_length*i
        handle.seek(offset)
        raw_data[zzz] = handle.read(row_length)
    raw_data = "".join(raw_data)

    # Parse the big string.
    all_format = ">%d%s" % (num_samples*len_indexes, typecode)
    all_data = struct.unpack(all_format, raw_data)

    # Split into matrix.
    data = [None] * len_indexes
    for i in range(len_indexes):
        j = i * num_samples
        data[i] = all_data[j:j+num_samples]
    return data


ALPHA_Z_CUTOFF = {}
def load_alpha_z_cutoff(geo_path, gseid, gplid):
    global ALPHA_Z_CUTOFF
    if not ALPHA_Z_CUTOFF:
        ALPHA_Z_CUTOFF = _load_alpha_z_cutoff_h(geo_path)
    ds = gseid, gplid
    assert ds in ALPHA_Z_CUTOFF, "missing z-cutoff: %s %s" % ds
    return ALPHA_Z_CUTOFF[ds]

def _load_alpha_z_cutoff_h(geo_path):
    import iolib
    import filelib

    filename = os.path.join(geo_path, "alpha_z_cutoff.txt")
    assert filelib.exists_nz(filename)
    ds2z = {}
    x = open(filename).read()
    for cols in iolib.split_tdf(x):
    #for d in filelib.read_row(filename, "gseid:s gplid:s z_cutoff:f"):
        GSEID, GPLID, z_cutoff = cols
        # Do not load the datasets where I could not find a good
        # cutoff.
        z_cutoff = float(z_cutoff)
        if z_cutoff < 1E-2:
            continue
        ds2z[(GSEID, GPLID)] = z_cutoff
        #ds2z[(d.gseid, d.gplid)] = d.z_cutoff
    return ds2z

def read_and_unpack(handle, format):
    import struct
    x = struct.unpack(format, handle.read(struct.calcsize(format)))
    return x

def pack_and_write(handle, format, *args):
    import struct
    x = struct.pack(format, *args)
    handle.write(x)

def read_and_unpack_str(handle):
    l, = read_and_unpack(handle, ">i")
    x, = read_and_unpack(handle, ">%ds" % l)
    return x

def pack_and_write_str(handle, s):
    pack_and_write(handle, ">i%ds" % len(s), len(s), s)



# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    import cgeolib
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in cgeolib.__dict__.keys():
        if not name.startswith("__"):
            this_module.__dict__[name] = cgeolib.__dict__[name]
