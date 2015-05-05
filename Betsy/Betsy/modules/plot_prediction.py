#plot_prediction.py
import os
from genomicode import mplgraph, filelib
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
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

    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_prediction_bar fails' % outfile
    )
    out_node = bie3.Data(rulebase.PredictionPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    loocv = ''
    if antecedents.data.attributes['loocv'] == 'yes':
        loocv = 'loocv'
    filename = ('prediction_' + original_file + '_' +
                antecedents.data.attributes['classify_alg'] + loocv + '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
