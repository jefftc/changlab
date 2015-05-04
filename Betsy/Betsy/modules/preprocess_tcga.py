#preprocess_tcga.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    out_attributes = get_out_attributes(out_attributes, in_data)
    TCGA_BIN = config.download_tcga
    command = ['python', TCGA_BIN, '--input', in_data.identifier, '--data',
               in_data.data.attributes['preprocess'], '--process_only',
               outfile]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    x = process.communicate()
    error_message = x[1]
    if error_message:
        raise ValueError(error_message)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for preprocess_tcga fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signalfile_tcga_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
