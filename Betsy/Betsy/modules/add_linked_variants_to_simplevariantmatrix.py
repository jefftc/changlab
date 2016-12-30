from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import math
        from genomicode import filelib
        from genomicode import jmath
        from genomicode import AnnotationMatrix
        from genomicode import SimpleVariantMatrix
        from Betsy import module_utils as mlib

        svm_node = in_data
        filelib.assert_exists_nz(svm_node.identifier)

        linked_file = mlib.get_user_option(
            user_options, "linked_variants_file", not_empty=True,
            check_file=True)
        
        # Read the variant file.
        SVM = SimpleVariantMatrix.read_as_am(svm_node.identifier)
        CHROM = SVM["______Chrom"]
        POS = SVM["______Pos"]
        POS = [int(x) for x in POS]
        all_coords = {}  # (chrom, pos) -> 1
        for x in zip(CHROM, POS):
            all_coords[x] = 1

        # Read the linked variant file.
        # Chrom  Pos  Perc Linked  p
        coord2info = {}  # (chrom, pos) -> d
        for d in filelib.read_row(linked_file, header=1):
            pos = int(d.Pos)
            if (d.Chrom, pos) not in all_coords:
                continue
            coord2info[(d.Chrom, pos)] = d

        # Align the linked annotations to the matrix.
        linked_headers = ["Perc Linked", "Score"]
        annotations = []
        for (chrom, pos) in zip(CHROM, POS):
            if (chrom, pos) not in coord2info:
                x = [""] * len(linked_headers)
                annotations.append(x)
                continue
            d = coord2info[(chrom, pos)]
            score = -10*math.log(float(d.p), 10)
            x = d.Perc_Linked, score
            assert len(x) == len(linked_headers)
            annotations.append(x)
        # Convert the headers and annotations to SVM format.
        linked_headers = ["Linkage______%s" % x for x in linked_headers]
        linked_annotations = jmath.transpose(annotations)

        # Make the new SimpleVariantMatrix.
        # Figure out where to put these annotations.
        INDEX = 4
        ## If Annovar exists, put after.
        #I = [i for (i, x) in enumerate(SVM.headers)
        #     if x.upper().startswith("ANNOVAR")]
        #if I:
        #    INDEX = max(INDEX, max(I)+1)
        headers = SVM.headers[:INDEX] + linked_headers + SVM.headers[INDEX:]
        x = [SVM.header2annots[x] for x in SVM.headers_h]
        all_annots = x[:INDEX] + linked_annotations + x[INDEX:]
        merged = AnnotationMatrix.create_from_annotations(
            headers, all_annots, headerlines=SVM.headerlines)

        SimpleVariantMatrix.write_from_am(outfile, merged)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
