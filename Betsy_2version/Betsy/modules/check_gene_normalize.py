#check_gene_normalize.py
import os
import shutil
from Betsy import module_utils,rulebase,bie
from genomicode import binreg
import arrayio
import numpy

def run(data_node,parameters, network):
    """check gene normalize"""
    outfile = name_outfile(data_node)
    parameters = get_out_attributes(parameters,data_node)
    shutil.copyfile(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for check_gene_normalize fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile1,**new_parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_check_normalize_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    M = arrayio.read(data_node.attributes['filename'])
    if is_gene_normalize_varaince(M):
        new_parameters['gene_normalize'] = 'variance'
    elif is_gene_normalize_ss(M):
        new_parameters['gene_normalize'] = 'sum_of_squares'
    else:
        new_parameters['gene_normalize'] = 'no'
    return new_parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def is_gene_normalize_varaince(M):
    for line in M.slice():
        if numpy.var(line)-1<0.000001:
            return False
    return True

def is_gene_normalize_ss(M):
    for line in M.slice():
        if numpy.sum([(x-numpy.mean(line))**2 for x in line])-1<0.000001:
            return False
    return True
