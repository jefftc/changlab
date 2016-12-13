"""Human readable format for variants.

Format:
                      Annovar          Coverage       <Sample>
                                                      <Caller>     <Caller>
Chrom  Pos  Ref  Alt  [Annovar Annots] <Samp> <Samp>  Ref/Alt/VAF  Ref/Alt/VAF


annot_matrix  Chrom, Pos, Annovar Annots, potentially other annotations.
named_matrix  Tuple of (<name>, AnnotationMatrix).  e.g., Coverage.
call_matrix   SparseCallMatrix.


Classes:
SimpleVariantMatrix
SparseCallMatrix
Call

Functions:
make_matrix
read
write

read_as_am      Return an AnnotationMatrix with samples and callers members.
write_from_am

"""
# _format_call
# _parse_call

class SimpleVariantMatrix:
    def __init__(self, samples, callers,
                 annot_matrix, named_matrices, call_matrix):
        # annot_matrix is an AnnotationMatrix.  It contains the
        # information for Chrom, Pos, Ref, Alt, and any Annovar
        # annotations.
        # call_matrix is a SparseCallMatrix that contains the Calls.
        from genomicode import AnnotationMatrix
        
        assert isinstance(annot_matrix, AnnotationMatrix.AnnotationMatrix)
        for name, matrix in named_matrices:
            assert isinstance(matrix, AnnotationMatrix.AnnotationMatrix)
            assert matrix.num_annots() == annot_matrix.num_annots()
        assert isinstance(call_matrix, SparseCallMatrix)

        # Figure out which column for sample and caller.
        ## assert len(matrix.headerlines) == 2
        ## header0 = matrix.headerlines[0].split("\t")
        ## header1 = matrix.headerlines[1].split("\t")
        ## assert len(header0) == len(header1)
        ## header_samples = [None] * len(header0)
        ## current_sample = None
        ## for i in range(1, len(header0)):
        ##     if header0[i]:
        ##         current_sample = header0[i]
        ##     header_samples[i] = current_sample
        ## header_callers = [None] * len(header1)
        ## for i in range(1, len(header1)):
        ##     if header1[i]:
        ##         header_callers[i] = header1[i]
        ## assert len(header_samples) == len(header_callers)
        ## samplecaller2i = {}  # (sample, caller) -> index
        ## for i, (s, c) in enumerate(zip(header_samples, header_callers)):
        ##     if s is None or c is None:
        ##         continue
        ##     assert s in samples, "Unknown sample: %s" % s
        ##     assert c in callers, "Unknown caller: %s" % c
        ##     assert (s, c) not in samplecaller2i
        ##     samplecaller2i[(s, c)] = i
        self.samples = samples[:]
        self.callers = callers[:]
        #self.annot_matrix = copy.deepcopy(annot_matrix)
        #self.call_matrix = copy.deepcopy(call_matrix)
        self.annot_matrix = annot_matrix
        self.named_matrices = named_matrices[:]
        self.call_matrix = call_matrix
    def num_variants(self):
        return self.annot_matrix.num_annots()


class SparseCallMatrix:
    def __init__(self, call_data):
        # call_data is a list of (chrom, pos, ref, alt, sample, caller, call).
        # chrom   string
        # pos     int
        # ref     string
        # alt     string
        # sample  string
        # caller  string
        # call    Call
        
        # (chrom, pos) -> (sample, caller) -> call
        coord2samplecaller2call = {}
        for x in call_data:
            chrom, pos, ref, alt, sample, caller, call = x
            assert type(pos) is type(0)
            assert isinstance(call, Call)
            coord = chrom, pos, ref, alt
            sc = sample, caller
            if coord not in coord2samplecaller2call:
                coord2samplecaller2call[coord] = {}
            # This can happen if the same sample was called in
            # multiple files (e.g. germline sample).
            x = "Possible that this sample was called in multiple files."
            assert sc not in coord2samplecaller2call[coord], \
                   "Duplicate: %s %s %s\n%s" % (caller, sample, coord, x)
            coord2samplecaller2call[coord][sc] = call
        self.coord2samplecaller2call = coord2samplecaller2call
    def has_call(self, chrom, pos, ref, alt, sample, caller):
        coord = chrom, pos, ref, alt
        sc = sample, caller
        if coord not in self.coord2samplecaller2call:
            return False
        return sc in self.coord2samplecaller2call[coord]
    def get_call(self, chrom, pos, ref, alt, sample, caller, default=None):
        # Return Call or default.
        coord = chrom, pos, ref, alt
        sc = sample, caller
        if coord not in self.coord2samplecaller2call:
            return default
        return self.coord2samplecaller2call[coord].get(sc, default)
    def set_call(self, chrom, pos, ref, alt, sample, caller, call):
        # call is typically Call.  To remove this call, set call to None.
        coord = chrom, pos, ref, alt
        sc = sample, caller

        if call is None:
            # Delete this call.
            if coord in self.coord2samplecaller2call:
                if sc in self.coord2samplecaller2call[coord]:
                    del self.coord2samplecaller2call[coord][sc]
                    if not self.coord2samplecaller2call[coord]:
                        del self.coord2samplecaller2call[coord]
            return
        if coord not in self.coord2samplecaller2call:
            self.coord2samplecaller2call[coord] = {}
        self.coord2samplecaller2call[coord][sc] = call


class Call:
    # num_ref   int or None
    # num_alt   int or None
    # vaf       float or None
    def __init__(self, num_ref, num_alt, vaf):
        assert num_ref is None or type(num_ref) is type(0)
        assert num_alt is None or type(num_alt) is type(0)
        assert vaf is None or type(vaf) is type(0.0)
        assert num_ref is not None or num_alt is not None or vaf is not None

        total = None
        if num_ref is not None and num_alt is not None:
            total = num_ref + num_alt
        self.num_ref = num_ref
        self.num_alt = num_alt
        self.vaf = vaf
        self.total = total


def make_matrix(samples, callers, annot_header, annot_data,
                named_data, call_data):
    # annot_header  list of headers for annot_data.
    # annot_data    list of tuples:  chrom, pos, ref, alt[, more]
    # named_data    list of (name, headers, all_annots)
    # call_data     list of tuples: chrom, pos, ref, alt, sample, caller, call
    # chrom   string
    # pos     int
    # ref     string
    # alt     string
    # sample  string
    # caller  string
    # call    Call object
    from genomicode import AnnotationMatrix

    # Make sure there's no duplicates.
    assert annot_header[:4] == ["Chrom", "Pos", "Ref", "Alt"]
    seen = {}
    for x in annot_data:
        x = x[:4]
        x = tuple(x)
        assert x not in seen, "Duplicate"
        seen[x] = 1

    # Make annotation matrix.
    for x in annot_data:
        assert len(x) == len(annot_header)
    headers = annot_header
    all_annots = []
    for i in range(len(headers)):
        x = [x[i] for x in annot_data]
        all_annots.append(x)
    annot_matrix = AnnotationMatrix.create_from_annotations(
        headers, all_annots)

    # Make named matrices.
    named_matrices = []
    for x in named_data:
        name, headers, all_annots = x
        matrix = AnnotationMatrix.create_from_annotations(
            headers, all_annots)
        x = name, matrix
        named_matrices.append(x)

    # Make call matrix.
    call_matrix = SparseCallMatrix(call_data)
    
    return SimpleVariantMatrix(
        samples, callers, annot_matrix, named_matrices, call_matrix)

    
def read(filename, is_csv=False):
    # Everything are strings.  No numeric conversion.
    from genomicode import filelib
    #from genomicode import AnnotationMatrix
    
    delimiter = "\t"
    if is_csv:
        delimiter = ","

    matrix = []
    for x in filelib.read_cols(filename, delimiter=delimiter):
        matrix.append(x)
        #if len(matrix) > 50000:  # DEBUG
        #    break
    assert len(matrix) >= 3      # at least 3 rows for the header
    for i in range(1, len(matrix)):
        assert len(matrix[i]) == len(matrix[0])
    assert len(matrix[0]) >= 4   # Chrom, Pos, Ref, Alt
    assert len(matrix[0]) >= 5, "No calls"
    
    header0 = matrix[0]
    header1 = matrix[1]
    header2 = matrix[2]
    #assert header0[0] == "Sample"
    #assert header1[0] == "Caller"
    assert header2[:4] == ["Chrom", "Pos", "Ref", "Alt"]

    # Make a list of all samples.
    I = [i for (i, x) in enumerate(header2) if x == "Ref/Alt/VAF"]
    assert I
    x = [header0[i] for i in I]
    #x = header0[1:]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    samples = x

    # Make a list of all callers.
    x = [header1[i] for i in I]
    #x = header1[1:]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    callers = x

    # Figure out where the annotations end.
    for i in range(1, len(header0)):
        if header0[i]:
            break
    else:
        raise AssertionError, "No calls"
    annot_end = i

    # Make the annotation matrix.
    annot_header = header2[:annot_end]
    annot_data = [x[:annot_end] for x in matrix[3:]]

    # Find the start coordinates of the named matrices.
    x = [i for (i, x) in enumerate(header0) if x]
    x = [i for i in x if i not in I]
    I_named = x  # list of start index of the named matrices.
    I_coord = []  # list of (start, end) of named matrices.
    for i in range(len(I_named)):
        i_start = I_named[i]
        if i+1 < len(I_named):
            i_end = I_named[i+1]
        else:
            i_end = I[0]
        I_coord.append((i_start, i_end))
    # Make the named matrices.
    named_data = []  # list of (name, named_header, named_annots)
    for (i_start, i_end) in I_coord:
        name = header0[i_start]
        assert name
        named_header = header2[i_start:i_end]
        M = [x[i_start:i_end] for x in matrix[3:]]
        named_annots = []
        for j in range(len(named_header)):
            x = [M[i][j] for i in range(len(M))]
            named_annots.append(x)
        x = name, named_header, named_annots
        named_data.append(x)

    # Make the call_data.
    call_data = []
    header_samples = [None] * len(header0)
    for i in I:
        if header0[i]:
            header_samples[i] = header0[i]
        else:
            header_samples[i] = header_samples[i-1]
        assert header_samples[i]
    header_callers = [None] * len(header1)
    for i in I:
        header_callers[i] = header1[i]
        assert header_callers[i]
    for i in range(3, len(matrix)):
        chrom, pos, ref, alt = matrix[i][:4]
        pos = int(pos)
        for j in I:
            sample, caller = header_samples[j], header_callers[j]
            if not matrix[i][j]:
                continue
            call = _parse_call(matrix[i][j])
            x = chrom, pos, ref, alt, sample, caller, call
            call_data.append(x)

    return make_matrix(
        samples, callers, annot_header, annot_data, named_data, call_data)
    

def write(handle_or_file, variant_matrix):
    import itertools
    
    handle = handle_or_file
    if type(handle) is type(""):
        handle = open(handle, 'w')

    annot_matrix = variant_matrix.annot_matrix
    call_matrix = variant_matrix.call_matrix

    # Make the headers for the annotations.
    header2 = annot_matrix.headers[:]
    #header0 = ["Sample"] + [""] * (len(header2)-1)
    #header1 = ["Caller"] + [""] * (len(header2)-1)
    header0 = [""] * len(header2)
    header1 = [""] * len(header2)
    # Make the headers for the named matrices.
    for (name, matrix) in variant_matrix.named_matrices:
        x0 = [name] + [""]*(len(matrix.headers)-1)
        x1 = [""]*len(matrix.headers)
        x2 = matrix.headers
        header0 += x0
        header1 += x1
        header2 += x2

    # Make the headers for the calls.
    for sample in variant_matrix.samples:
        #x2 = ["Num Ref", "Num Alt", "VAF"]
        #x1 = ["%s - %s" % (sample, caller)] + [""]*(len(x2)-1)
        x0 = [sample] + [""]*(len(variant_matrix.callers)-1)
        x1 = variant_matrix.callers
        x2 = ["Ref/Alt/VAF"]*len(variant_matrix.callers)
        header0 += x0
        header1 += x1
        header2 += x2

    # Print out the headers.
    assert len(header0) == len(header1)
    assert len(header0) == len(header2)
    print >>handle, "\t".join(header0)
    print >>handle, "\t".join(header1)
    print >>handle, "\t".join(header2)

    # Cache for convenience.
    num_calls = len(variant_matrix.samples) * len(variant_matrix.callers)
    sc2j = {}
    for j, (sample, caller) in enumerate(itertools.product(
        variant_matrix.samples, variant_matrix.callers)):
        sc2j[(sample, caller)] = j
    
    for i in range(annot_matrix.num_annots()):
        # Get the data from annot_matrix.
        row0 = []
        for h in annot_matrix.headers_h:
            row0.append(annot_matrix.header2annots[h][i])
        # Get the data from the named matrices.
        row1 = []
        for (name, matrix) in variant_matrix.named_matrices:
            for h in matrix.headers_h:
                row1.append(matrix.header2annots[h][i])
        # Get the data from the calls.
        row2 = [""] * num_calls
        chrom, pos, ref, alt = row0[:4]
        pos = int(pos)
        coord = chrom, pos, ref, alt
        sc2call = call_matrix.coord2samplecaller2call.get(coord, {})
        for sc, call in sc2call.iteritems():
            sample, caller = sc
            assert sample in variant_matrix.samples, sample
            assert caller in variant_matrix.callers, caller
            assert call
            j = sc2j[sc]
            row2[j] = _format_call(call)
        # Assemble the data.
        row = row0 + row1 + row2
        assert len(row) == len(header0)
        print >>handle, "\t".join(map(str, row))


def _format_call(call):
    # <ref>/<alt>/<vaf>
    # Some callers do not provide VAF.
    if not call:
        return ""
    ref = alt = vaf = ""
    if call.num_ref is not None:
        ref = str(call.num_ref)
    if call.num_alt is not None:
        alt = str(call.num_alt)
    if call.vaf is not None:
        vaf = "%.3f" % call.vaf
    return "%s/%s/%s" % (ref, alt, vaf)


def _parse_call(call_str):
    # Ref/Alt/VAF
    # Return Call object.
    if not call_str:
        return None
    x = call_str.split("/")
    assert len(x) == 3
    ref = alt = vaf = None
    if x[0]:
        ref = int(x[0])
    if x[1]:
        alt = int(x[1])
    if x[2]:
        vaf = float(x[2])
    return Call(ref, alt, vaf)


## def get_indexes_of_calls(var_matrix):
##     I = []
##     for i, h in enumerate(var_matrix.matrix.headers):
##         if h == "Ref/Alt/VAF":
##             I.append(i)
##     assert I
##     return I


def write_from_am(handle_or_file, svm_matrix):
    from genomicode import jmath

    headers0 = [None] * len(svm_matrix.headers)
    headers1 = [None] * len(svm_matrix.headers)
    headers2 = [None] * len(svm_matrix.headers)
    for i, header in enumerate(svm_matrix.headers):
        x = header.split("___")
        assert len(x) == 3, "Invalid header format: %s" % x
        headers0[i] = x[0]
        headers1[i] = x[1]
        headers2[i] = x[2]

    for i in range(len(headers0)-1, -1, -1):
        # If headers1[i] is the same as header1[i-1], then do not
        # write it out again.
        #
        # Exception: If headers0[i] != headers0[i-1], then we're
        # starting a new "block", and headers1[i] should still be
        # written out.
        # Example: If there's only one <Caller>, then the <Sample>
        # will not be blank, but the <Caller> should still be copied
        # over (because they are the same).
        # <Sample1>    <Sample2>
        # <Caller>     <Caller>
        # Ref/Alt/VAF  Ref/Alt/VAF
        if headers1[i] == headers1[i-1] and headers0[i] == headers0[i-1]:
            headers1[i] = ""
        if headers0[i] == headers0[i-1]:
            headers0[i] = ""

    matrix = []
    for i, header_h in enumerate(svm_matrix.headers_h):
        h0 = headers0[i]
        h1 = headers1[i]
        h2 = headers2[i]
        annots = svm_matrix.header2annots[header_h]
        x = [h0, h1, h2] + annots
        matrix.append(x)
    # Transpose the matrix.
    matrix = jmath.transpose(matrix)

    handle = handle_or_file
    if type(handle) is type(""):
        handle = open(handle, 'w')
    for x in svm_matrix.headerlines:
        print >>handle, x
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


def read_as_am(filename, is_csv=False):
    # Read file in SVM format.  Return an AnnotationMatrix object.
    # Does no special processing on any columns (i.e. no parsing as
    # integers or Call objects).  Everything is a string.
    
    # Header format:  <header0>___<header1>___<header2>
    # "blanks" are filled in.  E.g. "Annovar" occurs in each Annovar
    # column in header0.
    #
    # Headers:
    # ______Chrom
    # ______Pos
    # ______Ref
    # ______Alt
    # Num Callers______<Sample>
    # ...
    from genomicode import filelib
    from genomicode import AnnotationMatrix
    
    delimiter = "\t"
    if is_csv:
        delimiter = ","

    matrix = []
    for x in filelib.read_cols(filename, delimiter=delimiter):
        matrix.append(x)
    assert len(matrix) >= 3      # at least 3 rows for the header
    for i in range(1, len(matrix)):
        assert len(matrix[i]) == len(matrix[0])
    assert len(matrix[0]) >= 4   # Chrom, Pos, Ref, Alt
    assert len(matrix[0]) >= 5, "No calls"

    header0 = matrix[0]
    header1 = matrix[1]
    header2 = matrix[2]
    assert header2[:4] == ["Chrom", "Pos", "Ref", "Alt"]

    # Fill in the blanks for header1.
    for i in range(1, len(header1)):
        if header1[i]:
            continue
        # header1[i] is blank.  If header0[i], then this starts a new
        # "block".  Start with a new header1, and do not copy the old
        # one over.
        if not header1[i] and not header0[i]:
            header1[i] = header1[i-1]
    # Fill in the blanks for header0.
    for i in range(1, len(header0)):
        if not header0[i]:
            header0[i] = header0[i-1]

    # Make a list of all samples.
    I = [i for (i, x) in enumerate(header2) if x == "Ref/Alt/VAF"]
    assert I
    x = [header0[i] for i in I]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    samples = x

    # Make a list of all callers.
    x = [header1[i] for i in I]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    callers = x

    headers = []
    for x in zip(header0, header1, header2):
        x = "___".join(x)
        headers.append(x)
    all_annots = []
    for j in range(len(headers)):
        annots = [x[j] for x in matrix[3:]]
        all_annots.append(annots)
    matrix = AnnotationMatrix.create_from_annotations(headers, all_annots)
    matrix.samples = samples
    matrix.callers = callers
    return matrix
