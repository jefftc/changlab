from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from Betsy import module_utils
        in_data = antecedents
        assert module_utils.is_missing(in_data.identifier), 'no missing values'
        M = arrayio.read(in_data.identifier)
        f_out = file(outfile, 'w')
        for i in range(M.dim()[0]):
            for j in range(M.dim()[1]):
                if M._X[i][j] is None:
                    M._X[i][j] = '0'
    
        
        arrayio.tab_delimited_format.write(M, f_out)
        f_out.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for zero_fill_if_missing fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_zero_fill_' + original_file + '.tdf'
        return filename



