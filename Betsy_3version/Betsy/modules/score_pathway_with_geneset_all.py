#score_pathway_with_geneset_all.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(in_nodes, parameters, user_input, network):
    """analyze geneset"""
    data_node,geneset_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    score_geneset_path = config.score_geneset
    score_geneset_BIN = module_utils.which(score_geneset_path)
    assert score_geneset_BIN,'cannot find the %s' %score_geneset_path
    automatch = parameters['automatch']
    command = ['python', score_geneset_BIN, '-o', outfile,
               '--geneset_file', geneset_node.identifier,
               data_node.identifier,'--all']
    if automatch == 'yes':
        command.append('--automatch')
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for score_pathway_with_geneset fails' % outfile)
    out_node = bie3.Data(rulebase.GenesetAnalysis,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile')
    geneset_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='GenesetFile')
    return data_node, geneset_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'score_geneset_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)






