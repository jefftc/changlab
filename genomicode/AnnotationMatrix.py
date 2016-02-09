"""

Column-oriented tab-delimited text file.  First row is headers.  All
following rows contain data.

Classes:
AnnotationMatrix

Functions:
create_from_annot_matrix

uniquify_headers
replace_headers

colslice   Slice the columns based on a list of indexes.

read
write

"""

class AnnotationMatrix:
    def __init__(self, headers, headers_h, header2annots):
        # headers is a list of the original headers.
        # headers_h are the headers, hashed to ensure uniqueness.
        # headers2annots is a dictionary of hashed headers to the list
        # of annotations.
        assert headers
        assert headers_h
        assert len(headers) == len(headers_h)
        assert sorted(headers_h) == sorted(header2annots)
        for x in headers_h[1:]:
            assert len(header2annots[x]) == len(header2annots[headers_h[0]])
        self.headers = headers[:]
        self.headers_h = headers_h[:]
        self.header2annots = header2annots.copy()
    def get_annots(self, header):
        # Return a list of the annotations for this header.
        I = [i for i in range(len(self.headers)) if self.headers[i] == header]
        assert I, "header not found: %s" % header
        assert len(I) == 1, "Multiple headers match: %s" % header
        i = I[0]
        h = self.headers_h[i]
        return self.header2annots[h]
    def num_headers(self):
        return len(self.headers)
    def num_annots(self):
        if not self.headers_h:
            return 0
        h = self.headers_h[0]
        return len(self.header2annots[h])
    def copy(self):
        return AnnotationMatrix(
            self.headers, self.headers_h, self.header2annots)
    def normalize_header(self, header, index_base1=False):
        # Return the hashed header.  header may be either a header,
        # hashed header, or a 0-based index.  If index_base1 is True,
        # then indexes will be interpreted as 1-based.  If header
        # cannot be found, return None.
        if header in self.headers:
            i = self.headers.index(header)
            return self.headers_h[i]
        if header in self.headers_h:
            return header
        header_i = None
        try:
            header_i = int(header)
        except ValueError, x:
            pass
        if header_i is not None:
            if index_base1:
                header_i = head_i - 1
            assert header_i >= 0 and header_i < len(self.headers)
            return self.headers_h[i]
        return None
        #raise KeyError, header


def create_from_annotations(headers, all_annots):
    # headers is a list of headers.
    # all_annots is a list (parallel to headers) that contain the
    # annotations.
    assert headers
    assert len(headers) == len(all_annots)
    num_annots = len(all_annots[0])
    for x in all_annots[1:]:
        assert len(x) == num_annots
    
    headers_h = uniquify_headers(headers)
    assert len(headers_h) == len(all_annots)
    header2annots = {}
    for (header_h, annots) in zip(headers_h, all_annots):
        header2annots[header_h] = annots
    return AnnotationMatrix(headers, headers_h, header2annots)
    

def read(filename, is_csv=False, ignore_lines_startswith=None):
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
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True,
        delimiter=delimiter, ignore_lines_startswith=ignore_lines_startswith):
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
    return AnnotationMatrix(all_headers, headers_h, header2annots)


def write(handle_or_file, annot_matrix):
    from genomicode import jmath
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
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


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
    return AnnotationMatrix(new_headers, new_headers_h, header2annots)


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
    return AnnotationMatrix(headers, headers_h, header2annots)
