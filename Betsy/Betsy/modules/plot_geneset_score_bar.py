from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        import os
        from Betsy import module_utils
        from genomicode import jmath
        in_data = antecedents
        matrix = [x for x in filelib.read_cols(in_data.identifier)]
        matrix = [x[1:] for x in matrix]
        matrix = jmath.transpose(matrix)
        sample = matrix[0][1:]
        data = matrix[1:]
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        
        for one_data in data:
            value = one_data[1:]
            value = [float(i) for i in value]
            pair = [(value[i], sample[i]) for i in range(len(value))]
            pair.sort()
            gene_value = [i[0] for i in pair]
            label = [i[1] for i in pair]
            ylabel = one_data[0]
            from genomicode import mplgraph
            fig = mplgraph.barplot(gene_value,
                                   box_label=label,
                                   xtick_rotation=90,
                                   xlabel='sample',
                                   ylabel=ylabel)
            output = os.path.join(outfile, ylabel)
            fig.savefig(output + '.png')
    
        
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for plot_geneset_score_bar fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'geneset_plot_' + original_file + '.png'
        return filename



