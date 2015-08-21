from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """convert signal file to gct format"""
        import arrayio
        from Betsy import module_utils
        in_data = antecedents
        f = file(outfile, 'w')
        M = arrayio.read(in_data.identifier)
        M_c = arrayio.convert(M, to_format=arrayio.gct_format)
        arrayio.gct_format.write(M_c, f)
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for convert_signal_to_gct fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_' + original_file + '.gct'
        return filename



