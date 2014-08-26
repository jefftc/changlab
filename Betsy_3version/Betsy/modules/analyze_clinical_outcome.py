#analyze_clinical_outcome.py
import os
from Betsy import module_utils, bie3, rulebase
from genomicode import config
import subprocess
import shutil

def run(in_nodes, parameters, user_input, network):
    data_node, clinical_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    analyze_path = config.analyze_clinical
    analyze_BIN = module_utils.which(analyze_path)
    assert analyze_BIN,'cannot find the %s' %analyze_path
    x = []
    if 'rank_cutoff' in user_input:
        x.extend(['--rank_cutoff', user_input['rank_cutoff']])
    if 'zscore_cutoff' in user_input:
        x.extend(['--zscore_cutoff', user_input['zscore_cutoff']])
    command = ['python', analyze_BIN, data_node.identifier,
               clinical_node.identifier,'--outcome',
               user_input['outcome']+','+user_input['dead'],
               '--gene',user_input['genename'],'-o','clin'] + x
    process = subprocess.Popen(command, shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    outputfiles = os.listdir(os.getcwd())
    os.mkdir(outfile)
    for i in outputfiles:
        shutil.copyfile(i,os.path.join(outfile,i))
    assert module_utils.exists_nz(outfile),(
        'the output file %s for analyze_clinical_outcome fails'%outfile)
    out_node = bie3.Data(rulebase.ClinicalAnalysis,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

    
def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,clinical_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)


def name_outfile(in_nodes,user_input):
    data_node,clinical_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'survial_analysis_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile')
    clinical_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='ClinicalFile')
    return data_node, clinical_node
