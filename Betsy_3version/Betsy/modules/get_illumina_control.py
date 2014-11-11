#get_illumina_control.py

import shutil
import os
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, user_input,network,num_cores):
    outfile = name_outfile(data_node,user_input)
    result_files = os.listdir(data_node.identifier)
    for result_file in result_files:
        if '-controls' in result_file:
            goal_file = os.path.join(data_node.identifier,
                                     result_file)
            shutil.copyfile(goal_file,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for illu_control fails'%outfile)
    out_node = bie3.Data(rulebase.ControlFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'control_illumina_' + original_file +'.gct'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
