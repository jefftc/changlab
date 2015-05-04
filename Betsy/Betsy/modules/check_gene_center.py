#check_gene_center.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie3, rulebase
import numpy


def run(network, antecedents, out_attributes, user_options, num_cores):
    """check gene cetenr"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = get_out_attributes(out_attributes, in_data)
    shutil.copyfile(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for check_gene_center fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Normalize, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_check_center_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    new_parameters = out_attributes.copy()
    M = arrayio.read(antecedents.identifier)
    if is_gene_center_mean(M):
        new_parameters['gene_center'] = 'mean'
    elif is_gene_center_median(M):
        new_parameters['gene_center'] = 'median'
    else:
        new_parameters['gene_center'] = 'no'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def is_gene_center_mean(M):
    for line in M.slice():
        if numpy.mean(line) > 0.0000001:
            return False
    return 'mean'


def is_gene_center_median(M):
    for line in M.slice():
        if numpy.median(line) > 0.0000001:
            return False
    return 'median'
