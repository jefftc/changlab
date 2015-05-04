#cluster_genes.py
import os
import subprocess
from Betsy import module_utils, bie3, rulebase
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    """clustering the input file"""
    in_data = antecedents
    CLUSTER_BIN = config.cluster
    cluster = module_utils.which(CLUSTER_BIN)
    assert cluster, 'cannot find the %s' % CLUSTER_BIN
    distance_para = {'correlation': '1', 'euclidean': '7'}
    dist = distance_para[out_attributes['distance']]
    com_parameter = ["-g", dist, '-pg', '-e', '1']
    outfile = name_outfile(in_data, user_options)
    command = [CLUSTER_BIN, '-f', in_data.identifier, '-u', outfile]
    for i in com_parameter:
        command.append(i)
    
    process = subprocess.Popen(command,
                               shell=False,
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
        'the output file %s for cluster_genes_by_pca fails' % outfile
    )
    out_node = bie3.Data(rulebase.ClusterFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'cluster_file_' + original_file + '.cdt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
