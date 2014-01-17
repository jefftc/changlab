#check_gene_center.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie, rulebase
import numpy

def run(data_node,parameters, network):
    """check gene cetenr"""
    outfile = name_outfile(data_node)
    parameters = get_out_attributes(parameters,data_node)
    shutil.copyfile(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for check_gene_center fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile1,**new_parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_check_center_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    M = arrayio.read(data_node.attributes['filename'])
    if is_gene_center_mean(M):
        new_parameters['gene_center'] = 'mean'
    elif is_gene_center_median(M):
        new_parameters['gene_center'] = 'median'
    else:
        new_parameters['gene_center'] = 'no'
    return new_parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
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
