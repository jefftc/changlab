from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """select a num of genes according to num_features"""
        from Betsy import module_utils
        import arrayio
        in_data = antecedents
        f = file(outfile, 'w')
        num_features = 500
        if 'num_features_value' in user_options:
            num_features = int(user_options['num_features_value'])
    
        
        assert num_features > 0, 'the num_features should be >0'
        M = arrayio.read(in_data.identifier)
        M_c = M.matrix(range(0, num_features), None)
        arrayio.tab_delimited_format.write(M_c, f)
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for select_first_n_genes fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_select_n_' + original_file + '.tdf'
        return filename



