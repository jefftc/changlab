from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import mplgraph
        import arrayio
        from Betsy import module_utils
        from genomicode import jmath
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        data = jmath.transpose(M._X)
        tickname = M._col_names['_SAMPLE_NAME']
        fig = mplgraph.boxplot(data,
                               xlabel='Sample Name',
                               ylabel='Signal',
                               title='Signal Intensity',
                               box_label=tickname)
        fig.savefig(outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for plot_intensity_boxplot fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'intensity_' + original_file + '.png'
        return filename


