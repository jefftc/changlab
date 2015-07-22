#plot_sample_pca_with_predictions_wv.py
import os
from Betsy import module_utils
from Betsy import read_label_file
from genomicode import pcalib, genesetlib
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, classify_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    result_data = genesetlib.read_tdf(classify_node.identifier,
                                      preserve_spaces=True,
                                      allow_duplicates=True)
    for i in result_data:
        if i[0] == 'Predicted_class':
            legend = i[2]
    colors = ['r', 'b', 'g', 'y']
    legend_dict = {}
    for index, item in enumerate(legend):
        if item not in legend_dict:
            legend_dict[item] = [index]
        else:
            legend_dict[item].append(index)
    color = [''] * len(legend)
    for index, key in enumerate(legend_dict.keys()):
        c = colors[index]
        for i in legend_dict[key]:
            color[i] = c
    module_utils.plot_pca(data_node.identifier, outfile, color, legend)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_sample_pca_with_predictions fails' % outfile
    )
    out_node = bie3.Data(rulebase.PredictionPCAPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    data_node, classify_node = antecedents
    original_file = module_utils.get_inputid(classify_node.identifier)
    loocv = ''
    if classify_node.data.attributes['loocv'] == 'yes':
        loocv = 'loocv'
    filename = ('prediction_pca_plot' + original_file + '_' +
                classify_node.data.attributes['classify_alg'] + loocv + '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, classify_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='PcaAnalysis')
    filter2 = module_utils.AntecedentFilter(
        datatype_name='ClassLabelFile',
        classify_alg=out_attributes["classify_alg"])
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x
