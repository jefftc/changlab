from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import module_utils
        in_data = antecedents
        module_utils.plot_line_keywd(in_data.identifier, 'housekeeping', outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for plot_illu_housekeeping_line fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'housekeeping_plot_' + original_file + '.png'
        return filename
