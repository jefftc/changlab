from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import mplgraph
        from genomicode import filelib
        from Betsy import module_utils
        in_data = antecedents
        matrix = [x for x in filelib.read_cols(in_data.identifier)]
        header = matrix[0]
        index = header.index('Confidence')
        matrix = matrix[1:]
        confidence = [float(i[index]) for i in matrix]
        sample = [i[0] for i in matrix]
        if confidence == [''] * len(matrix) or 'Correct?' in header:
            index = header.index('Predicted_class')
            class_value = [i[index] for i in matrix]
            label_dict = dict()
            label_list = []
            i = -1
            for label in class_value:
                if label not in label_dict.keys():
                    i = i + 1
                    label_dict[label] = i
                label_list.append(label_dict[label])
            yticks = label_dict.keys()
            ytick_pos = [label_dict[i] for i in label_dict.keys()]
            fig = mplgraph.barplot(label_list,
                                   box_label=sample,
                                   ylim=(-0.5, 1.5),
                                   ytick_pos=ytick_pos,
                                   yticks=yticks,
                                   xtick_rotation='vertical',
                                   ylabel='Prediction',
                                   xlabel='Sample')
            fig.savefig(outfile)
        else:
            fig = mplgraph.barplot(confidence,
                                   box_label=sample,
                                   ylim=(-1.5, 1.5),
                                   xtick_rotation='vertical',
                                   ylabel='Prediction',
                                   xlabel='Sample')
            fig.savefig(outfile)

    
        
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for plot_prediction_bar fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        loocv = ''
        if antecedents.data.attributes['loocv'] == 'yes':
            loocv = 'loocv'
        
        filename = ('prediction_' + original_file + '_' +
                    antecedents.data.attributes['classify_alg'] + loocv + '.png')
        return filename



