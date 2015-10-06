"""

Column-oriented tab-delimited text file.  First row is headers.  All
following rows contain data.

Classes:
AnnotationMatrix

Functions:
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


def read(filename, is_csv=False):
    # Everything are strings.  No numeric conversion.
    import re
    from genomicode import genesetlib

    delimiter = "\t"
    if is_csv:
        delimiter = ","

    all_headers, all_annots = [], []
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True,
        delimiter=delimiter):
        name, description, annots = x

        # Hack: Some files contain special characters, which mess up
        # alignment. Fix this here.
        # na\xc3\xafve-WIBR3.5 hESC
        # na\xe2\x80\x9a\xc3\xa0\xc3\xb6\xe2\x88\x9a\xc3\xb2ve-C1.2 hiPSC
        annots = [re.sub("na\\W+ve", "naive", x) for x in annots]

        all_headers.append(name)
        all_annots.append(annots)

    headers_h = _hash_headers_unique(all_headers)
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


def _hash_headers_unique(headers):
    # Make sure the headers in all_headers is unique.
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
