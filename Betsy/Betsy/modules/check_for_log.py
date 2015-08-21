from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """log the input file"""
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        shutil.copyfile(in_data.identifier, outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for log_signal fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_log_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        import arrayio
        from genomicode import binreg
        new_parameters = out_attributes.copy()
        M = arrayio.read(antecedents.identifier)
        if binreg.is_logged_array_data(M):
            new_parameters['logged'] = 'yes'
        else:
            new_parameters['logged'] = 'no'
        
        return new_parameters



