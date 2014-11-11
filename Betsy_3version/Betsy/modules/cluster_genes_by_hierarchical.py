#cluster_genes.py
import os
import subprocess
from Betsy import module_utils
from Betsy import bie3, rulebase
from genomicode import config


def run(data_node,parameters, user_input, network,num_cores):
    """clustering the input file"""
    CLUSTER_BIN = config.cluster
    cluster_module = module_utils.which(CLUSTER_BIN)
    assert cluster_module, 'cannot find the %s' % CLUSTER_BIN
    distance_para = {'correlation': '1', 'euclidean': '7'}
    dist = distance_para[parameters['distance']]  
    com_parameter = ['-m', 's', '-e', '1', '-g', dist]
    outfile = name_outfile(data_node,user_input)
    command = [CLUSTER_BIN, '-f', data_node.identifier,
               '-u', outfile]
    for i in com_parameter:
        command.append(i)
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    result_files = os.listdir(os.getcwd())
    result_format = 'cdt'
    for result_file in result_files:
        if result_file.endswith(result_format):
            os.rename(result_file, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for cluster_genes_by_hierarchical fails' % outfile)
    out_node = bie3.Data(rulebase.ClusterFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'cluster_file_' + original_file + '.cdt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
