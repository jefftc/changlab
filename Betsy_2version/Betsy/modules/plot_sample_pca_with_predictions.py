#plot_sample_pca_with_predictions_wv.py
import os
from Betsy import module_utils
import shutil
from Betsy import read_label_file
from genomicode import pcalib,genesetlib
import arrayio
from Betsy import bie,rulebase

def run(in_nodes,parameters, network):
    data_node,classify_node = in_nodes
    outfile = name_outfile(in_nodes)
    result_data = genesetlib.read_tdf(classify_node.attributes['filename'],
                                      preserve_spaces=True,allow_duplicates=True)
    for i in result_data:
        if i[0] == 'Predicted_class':
            legend = i[2]
    colors = ['r','b','g','y']
    legend_dict = {}
    for index,item in enumerate(legend):
        if item not in legend_dict:
            legend_dict[item]=[index]
        else:
            legend_dict[item].append(index)
    color=['']*len(legend)
    for index,key in enumerate(legend_dict.keys()):
        c = colors[index]
        for i in legend_dict[key]:
            color[i]=c
    module_utils.plot_pca(data_node.attributes['filename'],outfile,color,legend)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_sample_pca_with_predictions fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.PredictionPCAPlot,**new_parameters)
    return out_node


def name_outfile(in_nodes):
    data_node,classify_node = in_nodes
    original_file = module_utils.get_inputid(classify_node.attributes['filename'])
    loocv = ''
    if classify_node.attributes['loocv'] == 'yes':
        loocv = 'loocv'
    filename = 'prediction_pca_plot' + original_file + '_' + classify_node.attributes['classify_alg'] +loocv+'.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,classify_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='PcaAnalysis')
    classify_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='ClassifyFile',
                                            optional_key='classify_alg',
                                            optional_value=parameters['classify_alg'])
    print classify_node
    return data_node,classify_node
