from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib

        simple_node, coverage_node = antecedents
        filelib.assert_exists_nz(simple_node.identifier)
        filelib.assert_exists_nz(coverage_node.identifier)

        ## Figure out if I'm adding coverage data from DNA or RNA.
        #in_attrs = simple_node.data.attributes
        #out_attrs = out_attributes
        #name = "with_rna_coverage"
        #assert name in in_attrs and name in out_attrs
        #is_rna_cov = False
        #if in_attrs[name] == "no" and out_attrs[name] == "yes":
        #    is_rna_cov = True

        add_coverage_to_svm(
            simple_node.identifier, coverage_node.identifier, out_filename,
            False)

        
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



def add_coverage_to_svm(svm_file, coverage_file, outfile, is_rna_cov):
    from genomicode import jmath
    from genomicode import filelib
    from genomicode import AnnotationMatrix
    from genomicode import SimpleVariantMatrix
    
    # Read the variant file.
    SVM = SimpleVariantMatrix.read(svm_file)
    AM = SVM.annot_matrix
    assert "Chrom" in AM.headers
    assert "Pos" in AM.headers
    CHROM = AM["Chrom"]
    POS = AM["Pos"]
    POS = [int(x) for x in POS]

    # Read the coverage matrix.
    # Chrom  Pos  <Sample>  [<Sample> ...]
    # Pos is 1-based.
    coord2sample2cov = {}  # (chrom, pos) -> sample -> ref/alt/vaf
    cov_samples = {}
    for d in filelib.read_row(coverage_file, header=1):
        coord = d.Chrom, int(d.Pos)
        if coord not in coord2sample2cov:
            coord2sample2cov[coord] = {}
        for i in range(2, len(d._header)):
            sample = d._header[i]
            cov = d._cols[i]
            if not cov:
                continue
            #coord2sample2cov[coord][sample] = int(cov)
            coord2sample2cov[coord][sample] = cov
            cov_samples[sample] = 1

    # Make sure the samples from the variant matrix can be found
    # in the coverage matrix.
    missing = [x for x in SVM.samples if x not in cov_samples]
    assert len(missing) < len(SVM.samples), (
        "SimpleVariantMatrix and coverage file have "
        "no common samples.")
    # If the samples aren't sequenced at high coverage, it's
    # possible they just don't have reads at these positions.  Be
    # a little lenient here, and accept the file if some of the
    # samples overlap.
    #x = missing
    #if len(x) > 5:
    #    x = x[:5] + ["..."]
    #msg = "Samples (%d) not found in coverage file: %s" % (
    #    len(missing), ", ".join(x))
    #assert not missing, msg
    # Report the coverage for the samples at the intersection.
    SAMPLES = [x for x in SVM.samples if x in cov_samples]

    # Align the matrix to the simple variant matrix.
    #matrix = [[None]*len(SVM.samples) for i in range(AM.num_annots())]
    matrix = [[None]*len(SAMPLES) for i in range(AM.num_annots())]
    for i in range(AM.num_annots()):
        coord = CHROM[i], POS[i]
        sample2cov = coord2sample2cov.get(coord, {})
        x = [sample2cov.get(x, "") for x in SAMPLES]
        #x = map(str, x)
        matrix[i] = x

    # Add the matrix back to the simple variant matrix.
    headers = SAMPLES
    all_annots = jmath.transpose(matrix)
    name = "Coverage"
    # If this is being used to add RNA coverage, use a different
    # name.
    if is_rna_cov:
        name = "RNA Coverage"
    x = AnnotationMatrix.create_from_annotations(headers, all_annots)
    SVM.named_matrices.append((name, x))

    # Write to file.
    SimpleVariantMatrix.write(outfile, SVM)
