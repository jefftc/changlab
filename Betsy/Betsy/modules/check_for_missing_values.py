from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        shutil.copyfile(in_data.identifier, outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_missing_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        from Betsy import module_utils
        new_parameters = out_attributes.copy()
        if module_utils.is_missing(antecedents.identifier):
            new_parameters['missing_values'] = 'yes'
        else:
            new_parameters['missing_values'] = 'no'
        
        return new_parameters



