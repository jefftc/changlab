#plot_sample_pca_wo_label.py
import os
import shutil
import arrayio
import matplotlib.cm as cm
from Betsy import bie
from Betsy import rulebase
from Betsy import read_label_file
from Betsy import module_utils

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    module_utils.plot_pca(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for pca_sample_plot fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.PcaPlot,**new_parameters)
    return out_node

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = ('Pca_' + original_file + '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='PcaAnalysis')

    return data_node




