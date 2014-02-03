#check_gene_center.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie3, rulebase
import numpy

def run(data_node,parameters, user_input, network):
    """check gene cetenr"""
    outfile = name_outfile(data_node,user_input)
    parameters = get_out_attributes(parameters,data_node)
    shutil.copyfile(data_node.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for check_gene_center fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile1,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_check_center_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    M = arrayio.read(data_node.identifier)
    if is_gene_center_mean(M):
        new_parameters['gene_center'] = 'mean'
    elif is_gene_center_median(M):
        new_parameters['gene_center'] = 'median'
    else:
        new_parameters['gene_center'] = 'no'
    return new_parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def is_gene_center_mean(M):
    for line in M.slice():
        if numpy.mean(line)>0.0000001:
            return False
    return 'mean'

def is_gene_center_median(M):
    for line in M.slice():
        if numpy.median(line)>0.0000001:
            return False
    return 'median'
