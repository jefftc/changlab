#cluster_genes.py
import os
import subprocess
from Betsy import module_utils, bie, rulebase

def run(data_node,parameters, network):
    """clustering the input file"""
    CLUSTER_BIN = 'cluster'
    distance_para = {'correlation': '1', 'euclidean': '7'}
    dist = distance_para[parameters['distance']]
    com_parameter = ["-g", dist, '-pg', '-e', '1']
    outfile = name_outfile(data_node)
    command = [CLUSTER_BIN, '-f', data_node.attributes['filename'],
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
    result_format = 'pca_gene.coords.txt'
    for result_file in result_files:
        if result_file.endswith(result_format):
            os.rename(result_file, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for cluster_genes_by_pca fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.ClusterFile,**new_parameters)
    return out_node


def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'cluster_file_' + original_file + '.cdt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
