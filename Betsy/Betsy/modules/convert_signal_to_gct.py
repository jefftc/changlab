from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """convert signal file to gct format"""
        import arrayio
        
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        M_c = arrayio.convert(M, to_format=arrayio.gct_format)
        
        f = file(outfile, 'w')
        arrayio.gct_format.write(M_c, f)
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_' + original_file + '.gct'
        return filename



