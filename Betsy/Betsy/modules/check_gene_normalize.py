#check_gene_normalize.py
import os
import shutil
from Betsy import module_utils, rulebase, bie3
from genomicode import binreg
import arrayio
import numpy


def run(network, antecedents, out_attributes, user_options, num_cores):
    """check gene normalize"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = get_out_attributes(out_attributes, in_data)
    shutil.copyfile(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for check_gene_normalize fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Normalize, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_check_normalize_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    new_parameters = out_attributes.copy()
    M = arrayio.read(antecedents.identifier)
    if is_gene_normalize_varaince(M):
        new_parameters['gene_normalize'] = 'variance'
    elif is_gene_normalize_ss(M):
        new_parameters['gene_normalize'] = 'sum_of_squares'
    else:
        new_parameters['gene_normalize'] = 'no'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
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
