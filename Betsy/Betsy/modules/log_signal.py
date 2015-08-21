from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """log the input file"""
        from Betsy import module_utils
        from genomicode import binreg
        import arrayio
        import math
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        assert not binreg.is_logged_array_data(M), 'the file is logged'
        for i in range(len(M._X)):
            for j in range(len(M._X[i])):
                if M._X[i][j] is not None:
                    if float(M._X[i][j]) < 1:
                        M._X[i][j] = 1
                    M._X[i][j] = math.log(float(M._X[i][j]), 2)
    
        
        f = file(outfile, 'w')
        M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
        arrayio.tab_delimited_format.write(M_c, f)
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for log_signal fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_log_' + original_file + '.tdf'
        return filename




