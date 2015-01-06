#score_pathway_with_scoresig.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    scoresig_path = config.scoresig
    scoresig_BIN = module_utils.which(scoresig_path)
    assert scoresig_BIN,'cannot find the %s' %scoresig_path
    command = ['python', scoresig_BIN, '-r', data_node.identifier, '-m', data_node.identifier, '-j', '20', '-o', outfile]
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for run_scoresig does not exists'% outfile)
    out_node = bie3.Data(rulebase.SignatureScore,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile')

    return data_node
    

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signature_score' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)
