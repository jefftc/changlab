#check_for_missing.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie
from Betsy import rulebase

def run(data_node, parameters, network):
    """log the input file"""
    parameters = get_out_attributes(parameters, data_node)
    outfile = name_outfile(data_node)
    shutil.copyfile(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for check_for_missing fails'%outfile)
    parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_missing_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    if module_utils.is_missing(data_node.attributes['filename']):
        new_parameters['missing_values'] = 'yes'
    else:
        new_parameters['missing_values'] = 'no'
    return new_parameters


def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

