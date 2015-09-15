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
        filename = 'signal_log_' + original_file + '.tdf'
        return filename


    def set_out_attributes(self, antecedents, out_attributes):
        import arrayio
        from genomicode import binreg

        attrs = out_attributes.copy()
        M = arrayio.read(antecedents.identifier)
        if binreg.is_logged_array_data(M):
            attrs['logged'] = 'yes'
        else:
            attrs['logged'] = 'no'
        return attrs



