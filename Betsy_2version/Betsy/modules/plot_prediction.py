#plot_prediction.py
import sys
import os
from genomicode import mplgraph,filelib
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    matrix=[x for x in filelib.read_cols(data_node.attributes['filename'])]
    header = matrix[0]
    index = header.index('Confidence')
    matrix=matrix[1:]
    confidence = [float(i[index]) for i in matrix]
    sample=[i[0] for i in matrix]
    if confidence==['']*len(matrix):
        index = header.index('Predicted_class')
        class_value = [i[index] for i in matrix]
        label_dict = dict()
        label_list = []
        i = -1
        for label in class_value:
            if  label not in label_dict.keys():
                i = i+1
                label_dict[label] = i
            label_list.append(label_dict[label])
        yticks=label_dict.keys()
        ytick_pos=[label_dict[i] for i in label_dict.keys()]
        fig=mplgraph.barplot(label_list,box_label=sample,
                   ylim=(-0.5,1.5),ytick_pos=ytick_pos,yticks=yticks,
                   xtick_rotation='vertical',ylabel='Prediction',xlabel='Sample')
        fig.savefig(outfile)
    else:
        fig=mplgraph.barplot(confidence,box_label=sample,
                   ylim=(-1.5,1.5),
                   xtick_rotation='vertical',ylabel='Prediction',xlabel='Sample')
        fig.savefig(outfile)
    
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_prediction_bar fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.PredictionPlot,**new_parameters)
    return out_node
 

    
def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    loocv = ''
    if data_node.attributes['loocv'] == 'yes':
        loocv = 'loocv'
    filename = 'prediction_' + original_file + '_' + data_node.attributes['classify_alg'] +loocv+'.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node


