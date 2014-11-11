#preprocess_tcga.py
import os
import shutil
from Betsy import module_utils
from Betsy import bie3, rulebase
import subprocess
from genomicode import config

def run(data_node,parameters, user_input,network,num_cores):
    outfile = name_outfile(data_node,user_input)
    parameters = get_out_attributes(parameters,data_node)
    TCGA_BIN = config.download_tcga
    command = ['python', TCGA_BIN,'--input',data_node.identifier,
               '--data', data_node.data.attributes['data'],
               '--process_only', outfile]
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    x = process.communicate()
    error_message = x[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for preprocess_tcga fails'%outfile)
    out_node = bie3.Data(rulebase._SignalFile_Postprocess,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signalfile_tcga_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,pool,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            pool,user_attributes)
    return data_node


