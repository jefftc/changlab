from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        from genomicode import filelib
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        shutil.copyfile(in_data.identifier, outfile)
        assert filelib.exists_nz(outfile), (
            'the output file %s for convert_postprocess_impute fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_file1_' + original_file + '.tdf'
        return filename



