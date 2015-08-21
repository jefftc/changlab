from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import module_utils
        import plot_sample_pca
        
        plot_sample_pca.plot_pca(antecedents.identifier, outfile)



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = ('Pca_' + original_file + '.png')
        return filename


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        # XXX does this need to be here?
        filter1 = module_utils.AntecedentFilter(datatype_name='PcaAnalysis')
        data_node = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1)
        return data_node
