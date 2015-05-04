#check_for_missing.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """log the input file"""
    in_data = antecedents
    out_attributes = get_out_attributes(out_attributes, in_data)
    outfile = name_outfile(in_data, user_options)
    shutil.copyfile(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for check_for_missing fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Impute, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_missing_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    new_parameters = out_attributes.copy()
    if module_utils.is_missing(antecedents.identifier):
        new_parameters['missing_values'] = 'yes'
    else:
        new_parameters['missing_values'] = 'no'
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
