#get_illumina_control.py

import shutil
import os
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    result_files = os.listdir(data_node.attributes['filename'])
    for result_file in result_files:
        if '-controls' in result_file:
            goal_file = os.path.join(data_node.attributes['filename'],
                                     result_file)
            shutil.copyfile(goal_file,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for illu_control fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'control_illumina_' + original_file +'.gct'
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
