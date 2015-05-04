#download_tcga.py
import os
from Betsy import module_utils
from genomicode import config
from Betsy import bie3, rulebase
import subprocess


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = get_out_attributes(out_attributes, in_data)
    TCGA_BIN = config.download_tcga
    assert 'disease' in user_options
    if 'date' in user_options:
        x = ['--date', user_options['date']]
    else:
        x = []
    
    command = ['python', TCGA_BIN, '--disease', user_options['disease'],
               '--data', out_attributes['preprocess'], '--download_only'] + x
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    
    result_files = os.listdir(os.getcwd())
    result_format = 'tar.gz'
    for result_file in result_files:
        if result_file.endswith(result_format):
            os.rename(result_file, outfile)

    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for download_tcga fails' % outfile
    )
    out_node = bie3.Data(rulebase.TCGAFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'tcga_file' + original_file + '.tar.gz'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    assert 'disease' in user_options
    x = ''
    if 'date' in user_options:
        x = '_' + user_options['date']
    identifier = user_options['disease'] + '_' + out_attributes['preprocess'] + x
    return identifier


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
