from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import plot_sample_pca
        
        plot_sample_pca.plot_pca(antecedents.identifier, outfile)



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = ('Pca_' + original_file + '.png')
        return filename


    
