"""Human readable format for variants.

Format:
Sample                                       <Sample>
Caller                                       <Caller>      <Caller>
Chrom  Pos  Ref  Alt  [Annovar Annotations]  Ref/Alt/VAF   Ref/Alt/VAF


Classes:
SimpleVariantMatrix
SparseCallMatrix
Call

Functions:
make_matrix
read
write

"""

class SimpleVariantMatrix:
    def __init__(self, samples, callers, annot_matrix, call_matrix):
        # annot_matrix is an AnnotationMatrix.  It contains the
        # information for Chrom, Pos, Ref, Alt, and any Annovar
        # annotations.
        # call_matrix is a SparseCallMatrix that contains the Calls.
        #import copy
        from genomicode import AnnotationMatrix
        
        assert isinstance(annot_matrix, AnnotationMatrix.AnnotationMatrix)
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
        self.call_matrix = call_matrix
        #self.matrix = copy.copy(matrix)
        #self.samplecaller2i = samplecaller2i
    def num_variants(self):
        return self.annot_matrix.num_annots()
    #def has_call(self, chrom, pos, sample, caller):
    #    return self.call_matrix.has_call(chrom, pos, sample, caller)
    #def get_call(self, chrom, pos, sample, caller):
    #    return self.call_matrix.get_call(chrom, pos, sample, caller)
    #def set_call(self, chrom, pos, sample, caller):
    #    return self.call_matrix.set_call(chrom, pos, sample, caller)


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


def make_matrix(samples, callers, annot_header, annot_data, call_data):
    # annot_header is a list of headers for annot_data.
    # annot_data is a list of tuples:  chrom, pos, ref, alt[, more]
    # call_data is a list of tuples: chrom, pos, ref, alt, sample, caller, call
    # chrom   string
    # pos     int
    # ref     string
    # alt     string
    # sample  string
    # caller  string
    # call    Call object
    #import itertools
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

    # Make call matrix.
    call_matrix = SparseCallMatrix(call_data)
    
    return SimpleVariantMatrix(samples, callers, annot_matrix, call_matrix)

    ## positions = {}  # (chrom, pos) -> (ref, alt)
    ## for x in coord_data:
    ##     chrom, pos, ref, alt = x
    ##     assert type(pos) is type(0)
    ##     positions[(chrom, pos)] = (ref, alt)

    ## position2data = {}   # (chrom, pos) -> list of (sample, caller, Call)
    ## for x in call_data:
    ##     chrom, pos, sample, caller, call = x
    ##     assert type(pos) is type(0)
    ##     key = chrom, pos
    ##     value = sample, caller, call
    ##     if key not in position2data:
    ##         position2data[key] = []
    ##     position2data[key].append(value)

    ## matrix = []
    ## header1 = ["Sample", "", "", ""]
    ## header2 = ["Caller", "", "", ""]
    ## header3 = ["Chrom", "Pos", "Ref", "Alt"]
    ## for sample in samples:
    ##     #x2 = ["Num Ref", "Num Alt", "VAF"]
    ##     #x1 = ["%s - %s" % (sample, caller)] + [""]*(len(x2)-1)
    ##     x1 = [sample] + [""]*(len(callers)-1)
    ##     x2 = callers
    ##     x3 = ["Ref/Alt/VAF"]*len(callers)
    ##     header1 += x1
    ##     header2 += x2
    ##     header3 += x3
    ## assert len(header1) == len(header2)
    ## assert len(header1) == len(header3)
    ## matrix.append(header1)
    ## matrix.append(header2)
    ## matrix.append(header3)

    ## for x in sorted(positions):
    ##     chrom, pos = x
    ##     ref, alt = positions[x]
    ##     assert x in position2data
    ##     data = position2data[x]

    ##     row = [chrom, pos, ref, alt]
    ##     for (sample, caller) in itertools.product(samples, callers):
    ##         # Look for the right object data lines.
    ##         for x in data:
    ##             s, c, call = x
    ##             if s == sample and c == caller:
    ##                 break
    ##         else:
    ##             call = None
    ##             row.append("")
    ##             continue
    ##         x = _format_call(call)
    ##         row.append(x)
    ##     assert len(row) == len(header1)
    ##     matrix.append(row)

    ## # Convert matrix to an AnnotationMatrix.
    ## headerlines = []
    ## headerlines.append("\t".join(matrix[0]))
    ## headerlines.append("\t".join(matrix[1]))
    ## headers = matrix[2]
    ## all_annots = []
    ## for i in range(len(headers)):
    ##     x = [x[i] for x in matrix[3:]]
    ##     all_annots.append(x)
    ## matrix = AnnotationMatrix.create_from_annotations(
    ##     headers, all_annots, headerlines)
    ## return SimpleVariantMatrix(samples, callers, matrix)

    
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
    assert len(matrix) >= 3
    for i in range(1, len(matrix)):
        assert len(matrix[i]) == len(matrix[0])
    assert len(matrix[0]) >= 5
    header0 = matrix[0]
    header1 = matrix[1]
    header2 = matrix[2]
    assert header0[0] == "Sample"
    assert header1[0] == "Caller"
    assert header2[:4] == ["Chrom", "Pos", "Ref", "Alt"]

    # Make a list of all samples.
    x = header0[1:]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    samples = x

    # Make a list of all callers.
    x = header1[1:]
    x = [x for x in x if x]
    # Get rid of duplicates, preserving order.
    x = [x[i] for (i, y) in enumerate(x) if y not in x[:i]]
    callers = x

    # Figure out where the annotations end and the calls start.
    for i in range(1, len(header0)):
        if header0[i]:
            break
    else:
        raise AssertionError, "No calls"
    call_start = i

    # Make the annotation matrix.
    annot_header = header2[:call_start]
    annot_data = [x[:call_start] for x in matrix[3:]]

    # Make the call_data.
    call_data = []
    header_samples = [None] * len(header0)
    for i in range(call_start, len(header0)):
        if header0[i]:
            header_samples[i] = header0[i]
        else:
            header_samples[i] = header_samples[i-1]
        assert header_samples[i]
    header_callers = [None] * len(header1)
    for i in range(call_start, len(header1)):
        header_callers[i] = header1[i]
        assert header_callers[i]
    for i in range(3, len(matrix)):
        chrom, pos, ref, alt = matrix[i][:4]
        pos = int(pos)
        for j in range(call_start, len(header0)):
            sample, caller = header_samples[j], header_callers[j]
            if not matrix[i][j]:
                continue
            call = _parse_call(matrix[i][j])
            x = chrom, pos, ref, alt, sample, caller, call
            call_data.append(x)

    return make_matrix(samples, callers, annot_header, annot_data, call_data)
    

def write(handle_or_file, variant_matrix):
    import itertools
    
    handle = handle_or_file
    if type(handle) is type(""):
        handle = open(handle, 'w')

    annot_matrix = variant_matrix.annot_matrix
    call_matrix = variant_matrix.call_matrix
    
    header2 = annot_matrix.headers[:]
    header0 = ["Sample"] + [""] * (len(header2)-1)
    header1 = ["Caller"] + [""] * (len(header2)-1)
    for sample in variant_matrix.samples:
        #x2 = ["Num Ref", "Num Alt", "VAF"]
        #x1 = ["%s - %s" % (sample, caller)] + [""]*(len(x2)-1)
        x0 = [sample] + [""]*(len(variant_matrix.callers)-1)
        x1 = variant_matrix.callers
        x2 = ["Ref/Alt/VAF"]*len(variant_matrix.callers)
        header0 += x0
        header1 += x1
        header2 += x2
    assert len(header0) == len(header1)
    assert len(header0) == len(header2)
    print >>handle, "\t".join(header0)
    print >>handle, "\t".join(header1)
    print >>handle, "\t".join(header2)

    # Cache for convenience.
    j2annots = {}
    for j, h in enumerate(annot_matrix.headers_h):
        annots = annot_matrix.header2annots[h]
        j2annots[j] = annots
    num_annots = len(j2annots)
    num_calls = len(variant_matrix.samples) * len(variant_matrix.callers)

    sc2j = {}
    for j, (sample, caller) in enumerate(itertools.product(
        variant_matrix.samples, variant_matrix.callers)):
        sc2j[(sample, caller)] = j
    
    for i in range(annot_matrix.num_annots()):
        row0 = [None] * num_annots
        for j in range(num_annots):
            row0[j] = j2annots[j][i]
        row1 = [""] * num_calls
        chrom, pos, ref, alt = row0[:4]
        pos = int(pos)
        coord = chrom, pos, ref, alt
        sc2call = call_matrix.coord2samplecaller2call.get(coord, {})
        for sc, call in sc2call.iteritems():
            assert call
            j = sc2j[sc]
            row1[j] = _format_call(call)
        row = row0 + row1
        assert len(row) == len(header0)
        print >>handle, "\t".join(map(str, row))


def _format_call(call):
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
    # Return tuple of ref, alt, vaf.  Each can be None if missing.
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
