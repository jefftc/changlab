#score_pathway_with_scoresig.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(in_nodes, parameters, user_input, network,num_cores):
    rma_node, mas5_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    scoresig_path = config.scoresig
    scoresig_BIN = module_utils.which(scoresig_path)
    assert scoresig_BIN,'cannot find the %s' %scoresig_path
    file1,file2 = module_utils.convert_to_same_platform(rma_node.identifier,
                                                        mas5_node.identifier)
    command = ['python', scoresig_BIN, '-r', file1, '-m', file2, '-j', '20', '-o', outfile]
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
    rma_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile',
                                           optional_key='preprocess',optional_value='rma')
    mas5_node = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='SignalFile',
                                            optional_key='preprocess',optional_value='mas5')
    return rma_node, mas5_node
    

def name_outfile(in_nodes,user_input):
    rma_node,mas5_node = in_nodes
    original_file = module_utils.get_inputid(
        rma_node.identifier)
    filename = 'signature_score' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    rma_node,mas5_node = in_nodes
    identifier = rma_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)
