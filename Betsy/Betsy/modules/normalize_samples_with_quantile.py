from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import quantnorm
        import arrayio
        from Betsy import module_utils
        from genomicode import filelib
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        Y = quantnorm.normalize(M)
        f = file(outfile, 'w')
        Y_c = arrayio.convert(Y, to_format=arrayio.pcl_format)
        arrayio.tab_delimited_format.write(Y_c, f)
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for quantile fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_quantile_' + original_file + '.tdf'
        return filename


    def set_out_attributes(self, antecedents, out_attributes):
        new_parameters = antecedents.data.attributes.copy()
        new_parameters['quantile_norm'] = 'yes'
        return new_parameters



