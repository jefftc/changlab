from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import AnnotationMatrix

        M = AnnotationMatrix.read(in_data.identifier, header_char="##")
        
        # Headers are:
        # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [Samples...]
        # Pull out the #CHROM and POS columns.
        assert M.num_headers()
        assert M.headers[0] == "#CHROM"
        assert M.headers[1] == "POS"
        chrom_annots = M["#CHROM"]
        pos_annots = M["POS"]

        lines = []
        seen = {}
        for chrom, pos in zip(chrom_annots, pos_annots):
            chrom, pos = chrom.strip(), pos.strip()
            x = chrom, pos
            if x in seen:
                continue
            seen[x] = 1
            x = "\t".join(x) + "\n"
            lines.append(x)
        open(out_filename, 'w').writelines(lines)


    def name_outfile(self, antecedents, user_options):
        return "positions.txt"
