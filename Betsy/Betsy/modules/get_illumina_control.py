#get_illumina_control.py

import shutil
import os
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    result_files = os.listdir(in_data.identifier)
    for result_file in result_files:
        if '-controls' in result_file:
            goal_file = os.path.join(in_data.identifier, result_file)
            shutil.copyfile(goal_file, outfile)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for illu_control fails' % outfile
    )
    out_node = bie3.Data(rulebase.ControlFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'control_illumina_' + original_file + '.gct'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
