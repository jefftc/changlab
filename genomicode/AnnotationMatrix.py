"""

Column-oriented tab-delimited text file.  First row is headers.  All
following rows contain data.

Classes:
AnnotationMatrix

Functions:
create_from_annotations
rowslice            Slice the rows based on a list of indexes.
colslice            Slice the columns based on a list of indexes.
replace_headers
uniquify_headers


read
write

"""
class AnnotationMatrix:
    def __init__(self, headers, headers_h, header2annots, headerlines=[]):
        # headers is a list of the original headers.
        # headers_h are the headers, hashed to ensure uniqueness.
        # headers2annots is a dictionary of hashed headers to the list
        # of annotations.  headerlines is a list of lines that occur
        # at the front of the file.  These might be the comment lines
        # that make up the "header" of VCF files.
        assert headers
        assert headers_h
        assert len(headers) == len(headers_h), "%d %d" % (
            len(headers), len(headers_h))
        assert sorted(headers_h) == sorted(header2annots)
        for x in headers_h[1:]:
            assert len(header2annots[x]) == len(header2annots[headers_h[0]])
        self.headers = headers[:]
        self.headers_h = headers_h[:]
        self.header2annots = header2annots.copy()
        for x in headerlines:
            assert type(x) is type("")
        self.headerlines = headerlines[:]  # no newlines
    #def get_annots(self, header):
    #    # Return a list of the annotations for this header.
    #    h = self.normalize_header(header)
    #    if annots is None:
    #        raise KeyError, header
    #    return self.header2annots[h]
    def num_headers(self):
        return len(self.headers)
    def num_annots(self):
        if not self.headers_h:
            return 0
        h = self.headers_h[0]
        return len(self.header2annots[h])
    def copy(self):
        return AnnotationMatrix(
            self.headers, self.headers_h, self.header2annots, self.headerlines)
    def __getitem__(self, header):
        # Return a list of the annotations for this header.
        h = self.normalize_header(header)
        if h is None:
            raise KeyError, header
        return self.header2annots[h]
    def __contains__(self, header):
        if self.normalize_header(header) is None:
            return False
        return True
    def normalize_header(self, header, index_base1=False):
        # Return the hashed header.  header may be either a header,
        # hashed header, or a 0-based index.  If index_base1 is True,
        # then indexes will be interpreted as 1-based.  If header
        # cannot be found, return None.

        # header can be:
        # 1.  Unique match to headers.
        # 2.  Unique match to headers_h.
        # 3.  Index.  (may be actual int or str)

        # Case 1.
        I = [i for i in range(len(self.headers)) if self.headers[i] == header]
        if len(I) == 1:
            return self.headers_h[I[0]]

        # Case 2.
        if header in self.headers_h:
            return header

        # Case 3.
        try:
            header_i = int(header)
        except ValueError, x:
            pass
        else:
            if index_base1:
                header_i = head_i - 1
            assert header_i >= 0 and header_i < len(self.headers)
            return self.headers_h[i]
        return None
        #raise KeyError, header
    def normalize_header_i(self, header, index_base1=False):
        # Return the 0-based index.  Same arguments as normalize_header.
        h = self.normalize_header(header, index_base1=index_base1)
        if h is None:
            return None
        return self.headers_h.index(h)


def create_from_annotations(headers, all_annots, headerlines=[]):
    # headers is a list of headers.
    # all_annots is a list (parallel to headers) that contain the
    # annotations.
    assert headers
    assert len(headers) == len(all_annots)
    num_annots = len(all_annots[0])
    for x in all_annots[1:]:
        assert len(x) == num_annots, "Inconsistent annot lengths"
    
    headers_h = uniquify_headers(headers)
    assert len(headers_h) == len(all_annots)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    x = AnnotationMatrix(
        headers, headers_h, header2annots, headerlines=headerlines)
    return x
    

def colslice(MATRIX, I):
    for i in I:
        assert i >= 0 and i < len(MATRIX.headers_h)
    new_headers = [MATRIX.headers[i] for i in I]
    old_headers_h = [MATRIX.headers_h[i] for i in I]
    new_headers_h = uniquify_headers(new_headers)
    header2annots = {}
    for i in range(len(old_headers_h)):
        oh, nh = old_headers_h[i], new_headers_h[i]
        header2annots[nh] = MATRIX.header2annots[oh]
    x = AnnotationMatrix(
        new_headers, new_headers_h, header2annots, MATRIX.headerlines)
    return x


def rowslice(MATRIX, I):
    num_annots = MATRIX.num_annots()
    for i in I:
        assert i >= 0 and i < num_annots
    header2annots = {}
    for header, old_annots in MATRIX.header2annots.iteritems():
        new_annots = [old_annots[i] for i in I]
        header2annots[header] = new_annots
    x = AnnotationMatrix(
        MATRIX.headers, MATRIX.headers_h, header2annots, MATRIX.headerlines)
    return x


def replace_headers(MATRIX, headers):
    # Return a new AnnotationMatrix with these headers.
    assert len(headers) == len(MATRIX.headers)
    headers_h = uniquify_headers(headers)
    header2annots = {}
    for header_old in MATRIX.header2annots:
        # Use the index to get the hashed header.
        i = MATRIX.headers_h.index(header_old)
        header_new = headers_h[i]
        header2annots[header_new] = MATRIX.header2annots[header_old]
    x = AnnotationMatrix(headers, headers_h, header2annots, MATRIX.headerlines)
    return x


def uniquify_headers(headers):
    # Return a parallel list of unique headers.
    header2I = {}  # header -> list of indexes
    for i, header in enumerate(headers):
        if header not in header2I:
            header2I[header] = []
        header2I[header].append(i)

    nodup = headers[:]
    for (header, I) in header2I.iteritems():
        if len(I) < 2:
            continue
        for i in range(len(I)):
            nodup[I[i]] = "%s_%d" % (header, i+1)
    return nodup


def read(filename, is_csv=False, header_char=None):
    # Everything are strings.  No numeric conversion.
    import re
    from genomicode import genesetlib

    delimiter = "\t"
    if is_csv:
        delimiter = ","

    # re.sub takes a lot of time (25% of all running time!).  Compile
    # it.
    re_naive = re.compile("na\\W+ve")

    all_headers, all_annots = [], []
    all_comments = []
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True,
        delimiter=delimiter, yield_lines_startswith=header_char):
        if type(x) is type(""):
            all_comments.append(x)
            continue
        name, description, annots = x

        # Hack: Some files contain special characters, which mess up
        # alignment. Fix this here.
        # na\xc3\xafve-WIBR3.5 hESC
        # na\xe2\x80\x9a\xc3\xa0\xc3\xb6\xe2\x88\x9a\xc3\xb2ve-C1.2 hiPSC
        #annots = [re.sub("na\\W+ve", "naive", x) for x in annots]
        annots = [re_naive.sub("naive", x) for x in annots]

        all_headers.append(name)
        all_annots.append(annots)
    assert all_headers, "Empty file: %s" % filename

    headers_h = uniquify_headers(all_headers)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    return AnnotationMatrix(
        all_headers, headers_h, header2annots, headerlines=all_comments)


def write(handle_or_file, annot_matrix, delim=None):
    from genomicode import jmath

    if delim is None:
        delim = "\t"
    
    matrix = []
    for i, header_h in enumerate(annot_matrix.headers_h):
        header = annot_matrix.headers[i]
        annots = annot_matrix.header2annots[header_h]
        x = [header] + annots
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)

    handle = handle_or_file
    if type(handle) is type(""):
        handle = open(handle, 'w')
    for x in annot_matrix.headerlines:
        print >>handle, x
    for x in matrix:
        print >>handle, delim.join(map(str, x))
