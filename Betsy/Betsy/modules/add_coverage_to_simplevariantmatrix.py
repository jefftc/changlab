from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        from genomicode import AnnotationMatrix
        from genomicode import SimpleVariantMatrix

        simple_node, coverage_node = antecedents
        filelib.assert_exists_nz(simple_node.identifier)
        filelib.assert_exists_nz(coverage_node.identifier)

        # Read the variant file.
        SVM = SimpleVariantMatrix.read(simple_node.identifier)
        AM = SVM.annot_matrix
        assert "Chrom" in AM.headers
        assert "Pos" in AM.headers
        CHROM = AM["Chrom"]
        POS = AM["Pos"]
        POS = [int(x) for x in POS]

        # Read the coverage matrix.
        # Chrom  Pos  <Sample>  [<Sample> ...]
        # Pos is 1-based.
        coord2sample2cov = {}  # (chrom, pos) -> sample -> coverage
        cov_samples = {}
        for d in filelib.read_row(coverage_node.identifier, header=1):
            coord = d.Chrom, int(d.Pos)
            if coord not in coord2sample2cov:
                coord2sample2cov[coord] = {}
            for i in range(2, len(d._header)):
                sample = d._header[i]
                cov = d._cols[i]
                if not cov:
                    continue
                coord2sample2cov[coord][sample] = int(cov)
                cov_samples[sample] = 1

        # Make sure the samples from the variant matrix can be found
        # in the coverage matrix.
        missing = [x for x in SVM.samples if x not in cov_samples]
        assert len(missing) < len(SVM.samples), (
            "SimpleVariantMatrix and coverage file have "
            "no common samples.")
        # If the samples aren't sequenced at high coverage, it's
        # possible they just don't have reads at these positions.  Be
        # a little lenient here, and accept the file is some of the
        # samples overlap.
        #x = missing
        #if len(x) > 5:
        #    x = x[:5] + ["..."]
        #msg = "Samples (%d) not found in coverage file: %s" % (
        #    len(missing), ", ".join(x))
        #assert not missing, msg

        # Align the matrix to the simple variant matrix.
        matrix = [[None]*len(SVM.samples) for i in range(AM.num_annots())]
        for i in range(AM.num_annots()):
            coord = CHROM[i], POS[i]
            sample2cov = coord2sample2cov.get(coord, {})
            x = [sample2cov.get(x, "") for x in SVM.samples]
            x = map(str, x)
            matrix[i] = x

        # Add the matrix back to the simple variant matrix.
        headers = SVM.samples
        all_annots = []
        for j in range(len(headers)):
            x = [matrix[i][j] for i in range(len(matrix))]
            all_annots.append(x)
        x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        SVM.named_matrices.append(("Coverage", x))

        # Write to file.
        SimpleVariantMatrix.write(out_filename, SVM)

        
        #x = ["%s Cov" % x for x in SVM.samples]
        #headers = AM.headers + x
        #all_annots = []
        #for h in AM.headers_h:
        #    all_annots.append(AM.header2annots[h])
        #for i in range(len(SVM.samples)):
        #    x = [x[i] for x in matrix]
        #    assert len(x) == AM.num_annots()
        #    all_annots.append(x)
        #assert not AM.headerlines
        #x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        #x = SimpleVariantMatrix.SimpleVariantMatrix(
        #    SVM.samples, SVM.callers, x, SVM.call_matrix)
        #SimpleVariantMatrix.write(out_filename, x)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
