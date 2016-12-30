from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from genomicode import jmath
        from genomicode import AnnotationMatrix
        from genomicode import SimpleVariantMatrix
        from Betsy import module_utils as mlib

        svm_node = in_data
        filelib.assert_exists_nz(svm_node.identifier)

        cosmic_file = mlib.get_user_option(
            user_options, "cosmic_variants_file", not_empty=True,
            check_file=True)
        
        # Read the variant file.
        SVM = SimpleVariantMatrix.read_as_am(svm_node.identifier)
        CHROM = SVM["______Chrom"]
        POS = SVM["______Pos"]
        POS = [int(x) for x in POS]
        all_coords = {}  # (chrom, pos) -> 1
        for x in zip(CHROM, POS):
            all_coords[x] = 1

        # Read the COSMIC variant file.
        # Chrom  Start  End  GRCh  Count  SNP
        # Mutation CDS  Mutation AA
        # FATHMM prediction  FATHMM score  Mutation somatic status
        coord2info = {}  # (chrom, pos) -> d
        for d in filelib.read_row(cosmic_file, header=1):
            start, end = int(d.Start), int(d.End)
            in_svm = False
            for pos in range(start, end+1):
                if (d.Chrom, pos) in all_coords:
                    in_svm = True
                    break
            if not in_svm:
                continue
            coord2info[(d.Chrom, pos)] = d

        # Align the COSMIC annotations to the matrix.
        cosmic_headers = [
            "SNP", "Num Tumors", "Mutation CDS", "Mutation AA",
            "FATHMM prediction", "FATHMM score", "Mutation somatic status"]
        annotations = []
        for (chrom, pos) in zip(CHROM, POS):
            if (chrom, pos) not in coord2info:
                x = [""] * len(cosmic_headers)
                annotations.append(x)
                continue
            d = coord2info[(chrom, pos)]
            x = d.SNP, d.Count, d.Mutation_CDS, d.Mutation_AA, \
                d.FATHMM_prediction, d.FATHMM_score, \
                d.Mutation_somatic_status
            annotations.append(x)
        # Convert the headers and annotations to SVM format.
        cosmic_headers = ["COSMIC______%s" % x for x in cosmic_headers]
        cosmic_annotations = jmath.transpose(annotations)

        # Make the new SimpleVariantMatrix.
        # Figure out where to put these annotations.
        INDEX = 4
        # If Annovar exists, put after.
        I = [i for (i, x) in enumerate(SVM.headers)
             if x.upper().startswith("ANNOVAR")]
        if I:
            INDEX = max(INDEX, max(I)+1)
        # If SnpEff exists, put after.
        I = [i for (i, x) in enumerate(SVM.headers)
             if x.upper().startswith("SNPEFF")]
        if I:
            INDEX = max(INDEX, max(I)+1)
        headers = SVM.headers[:INDEX] + cosmic_headers + SVM.headers[INDEX:]
        x = [SVM.header2annots[x] for x in SVM.headers_h]
        all_annots = x[:INDEX] + cosmic_annotations + x[INDEX:]
        merged = AnnotationMatrix.create_from_annotations(
            headers, all_annots, headerlines=SVM.headerlines)

        SimpleVariantMatrix.write_from_am(outfile, merged)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
