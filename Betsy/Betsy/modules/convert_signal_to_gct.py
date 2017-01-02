from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import arrayplatformlib as apl

        MATRIX = arrayio.read(antecedents.identifier)

        # Converter will just use the first two columns for NAME and
        # Description.  Try to find better ones.
        # Column 1: PROBE_ID, GENE_ID, <col 1>
        # Column 2: GENE_SYMBOL, <col 2>
        num_headers = len(MATRIX.row_names())
        cat2header = apl.categorize_headers(MATRIX, remove_version=True)

        h1 = cat2header.get(apl.PROBE_ID)
        if not h1:
            h1 = cat2header.get(apl.GENE_ID)
        if not h1 and num_headers:
            h1 = MATRIX.row_names()[0]

        h2 = cat2header.get(apl.GENE_SYMBOL)
        if not h2 and num_headers >= 2:
            h2 = MATRIX.row_names()[1]

        if h2:
            i = MATRIX._row_order.index(h2)
            MATRIX._row_order.pop(i)
            MATRIX._row_order.insert(0, h2)
        if h1:
            i = MATRIX._row_order.index(h1)
            MATRIX._row_order.pop(i)
            MATRIX._row_order.insert(0, h1)
        
        MATRIX_c = arrayio.convert(MATRIX, to_format=arrayio.gct_format)
        
        arrayio.gct_format.write(MATRIX_c, open(outfile, 'w'))


    def name_outfile(self, antecedents, user_options):
        return "signal.gct"



