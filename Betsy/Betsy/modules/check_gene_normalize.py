#check_gene_normalize.py
import os
import shutil
from Betsy import module_utils, rulebase, bie3
from genomicode import binreg
import arrayio
import numpy


def run(data_node, parameters, user_input, network, num_cores):
    """check gene normalize"""
    outfile = name_outfile(data_node, user_input)
    parameters = get_out_attributes(parameters, data_node)
    shutil.copyfile(data_node.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for check_gene_normalize fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Normalize, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_check_normalize_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_node):
    new_parameters = parameters.copy()
    M = arrayio.read(data_node.identifier)
    if is_gene_normalize_varaince(M):
        new_parameters['gene_normalize'] = 'variance'
    elif is_gene_normalize_ss(M):
        new_parameters['gene_normalize'] = 'sum_of_squares'
    else:
        new_parameters['gene_normalize'] = 'no'
    return new_parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node


def is_gene_normalize_varaince(M):
    for line in M.slice():
        if abs(numpy.var(line) - 1) > 0.000001:
            return False
    return True


def is_gene_normalize_ss(M):
    for line in M.slice():
        if abs(numpy.sum([(x - numpy.mean(line)) ** 2
                          for x in line]) - 1) > 0.000001:
            return False
    return True
