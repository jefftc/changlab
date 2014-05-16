#check_for_missing.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie3
from Betsy import rulebase

def run(data_node, parameters, user_input, network):
    """log the input file"""
    parameters = get_out_attributes(parameters, data_node)
    outfile = name_outfile(data_node,user_input)
    shutil.copyfile(data_node.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for check_for_missing fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile_Impute,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_missing_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    if module_utils.is_missing(data_node.identifier):
        new_parameters['missing_values'] = 'yes'
    else:
        new_parameters['missing_values'] = 'no'
    return new_parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

