from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        from Betsy import module_utils
        from genomicode import filelib
    
        in_data = antecedents

        # Copy the file objects so that gzip'd files get properly uncompressed.
        fsrc = filelib.openfh(in_data.identifier)
        fdst = open(outfile, 'w')
        shutil.copyfileobj(fsrc, fdst)
        fdst.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        format_type = antecedents.data.attributes['format']
        filename = 'signal_file_' + original_file + '.' + format_type
        return filename



