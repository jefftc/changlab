from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import SimpleVariantMatrix

        matrix = SimpleVariantMatrix.read(in_data.identifier)
        AM = matrix.annot_matrix
        assert "Chrom" in AM.headers
        assert "Pos" in AM.headers
        CHROM = AM["Chrom"]
        POS = AM["Pos"]
        POS = [int(x) for x in POS]

        # Write a positions file.
        # <chrom> <pos 0-based>
        #
        # SimpleVariantMatrix uses 1-based coordinates.  Change to
        # 0-based.

        handle = open(out_filename, 'w')
        seen = {}
        for i in range(AM.num_annots()):
            chrom = CHROM[i]
            pos = POS[i] - 1
            if (chrom, pos) in seen:
                continue
            seen[(chrom, pos)] = 1
            x = [chrom, pos]
            print >>handle, "\t".join(map(str, x))
        handle.close()
            

    def name_outfile(self, antecedents, user_options):
        return "positions.txt"
