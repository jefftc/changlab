#convert_postprocess_impute.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = set_out_attributes(in_data, out_attributes)
    shutil.copyfile(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_postprocess_impute fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Impute, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_file1_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
