#check_for_log.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """log the input file"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = get_out_attributes(out_attributes, in_data)
    shutil.copyfile(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for log_signal fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_log_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    new_parameters = out_attributes.copy()
    M = arrayio.read(antecedents.identifier)
    if binreg.is_logged_array_data(M):
        new_parameters['logged'] = 'yes'
    else:
        new_parameters['logged'] = 'no'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
