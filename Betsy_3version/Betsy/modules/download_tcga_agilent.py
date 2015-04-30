#download_tcga_agilent.py
import os
import shutil
from Betsy import module_utils
from genomicode import config
import arrayio
from Betsy import bie3, rulebase
import subprocess


def run(data_node, parameters, user_input, network, num_cores):
    outfile = name_outfile(data_node, user_input)
    parameters = get_out_attributes(parameters, data_node)
    TCGA_BIN = config.download_tcga
    assert 'disease' in user_input
    if 'date' in user_input:
        x = ['--date', user_input['date']]
    else:
        x = []
    command = ['python', TCGA_BIN, '--disease', user_input['disease'],
               '--data', parameters['preprocess'], '--download_and_extract'] + x
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    result_files = os.listdir(os.getcwd())
    result_format = '.txt'
    for result_file in result_files:
        if result_file.endswith(result_format) and result_file != 'MANIFEST.txt':
            os.rename(result_file, outfile)

    assert module_utils.exists_nz(outfile), (
        'the output file %s for download_tcga_agilent fails' % outfile
    )
    out_node = bie3.Data(rulebase.TCGAFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'tcga_file' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_node):
    return parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    assert 'disease' in user_input
    x = ''
    if 'date' in user_input:
        x = '_' + user_input['date']
    identifier = user_input['disease'] + '_' + parameters['preprocess'] + x
    return identifier


def find_antecedents(network, module_id, pool, parameters, user_attributes):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
