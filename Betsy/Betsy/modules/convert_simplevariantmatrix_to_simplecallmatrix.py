from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        #from genomicode import filelib
        from genomicode import SimpleVariantMatrix
        from Betsy import module_utils as mlib

        simple_file = in_data.identifier
        metadata = {}

        num_callers = mlib.get_user_option(
            user_options, "num_callers", not_empty=True, type=int)
        assert num_callers >= 0 and num_callers < 100

        var_matrix = SimpleVariantMatrix.read(simple_file)
        annot_matrix = var_matrix.annot_matrix
        call_matrix = var_matrix.call_matrix

        # For each coord and sample, count the number of callers.
        coord2sample2nc = {}  # (chrom, pos, ref, alt) -> sample -> num callers
        for x in call_matrix.coord2samplecaller2call.iteritems():
            coord, samplecaller2call = x
            if coord not in coord2sample2nc:
                coord2sample2nc[coord] = {}
            sample2nc = coord2sample2nc[coord]
            for (sample, caller), call in samplecaller2call.iteritems():
                # Make sure this is a real call.
                if not (call.num_ref or call.num_alt or
                        call.total or call.vaf):
                    continue
                sample2nc[sample] = sample2nc.get(sample, 0) + 1

        # Make a list of the coordinates that have the right number of calls.
        calls = {}  # coord -> sample -> nc
        for coord, sample2nc in coord2sample2nc.iteritems():
            for sample, nc in sample2nc.iteritems():
                if nc < num_callers:
                    continue
                if coord not in calls:
                    calls[coord] = {}
                calls[coord][sample] = nc

        handle = open(out_filename, 'w')

        # Print out the matrix.
        header = annot_matrix.headers + var_matrix.samples
        print >>handle, "\t".join(header)
        
        # Cache for convenience.
        j2annots = {}
        for j, h in enumerate(annot_matrix.headers_h):
            annots = annot_matrix.header2annots[h]
            j2annots[j] = annots
        num_annots = len(j2annots)

        chrom, pos = annot_matrix["Chrom"], annot_matrix["Pos"]
        ref, alt = annot_matrix["Ref"], annot_matrix["Alt"]
        pos = [int(x) for x in pos]
        for i, coord in enumerate(zip(chrom, pos, ref, alt)):
            if coord not in calls:
                continue

            row0 = [None] * num_annots
            for j in range(num_annots):
                row0[j] = j2annots[j][i]
            row1 = [""] * len(var_matrix.samples)
            for j, sample in enumerate(var_matrix.samples):
                if sample in calls[coord]:
                    row1[j] = coord2sample2nc[coord][sample]

            row = row0 + row1
            assert len(row) == len(header)
            print >>handle, "\t".join(map(str, row))


        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "calls.txt"
