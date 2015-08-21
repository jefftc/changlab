from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """unlog the pcl file"""
        import arrayio
        from genomicode import binreg

        M = arrayio.read(antecedents.identifier)
        assert binreg.is_logged_array_data(M), (
            'the input file %s should be logged' % antecedents.identitifer)
        
        for i in range(len(M._X)):
            for j in range(len(M._X[i])):
                if M._X[i][j] is not None:
                    M._X[i][j] = 2 ** float(M._X[i][j])
    
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M, f)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'signal_unlog' + original_file + '.tdf'
