#consensusClustering.py
import os
from Betsy import bie3,rule_engine_bie3
from Betsy import rulebase
from Betsy import module_utils
from genomicode import config
import subprocess


def run(data_node,parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    module_name = 'ConsensusClustering'
    gp_parameters = dict()
    file1 = data_node.identifier
    
    gp_parameters['input.filename'] = file1
    
    if 'cc_kmax' in user_input:
        assert module_utils.is_number(
            user_input['cc_kmax']), 'the cc_kmax should be number'
        gp_parameters['cc_kmax'] = str(user_input['cc_kmax'])
    if 'cc_resampling_iter' in user_input:
        assert module_utils.is_number(
            user_input['cc_resampling_iter']), 'the cc_resampling_iter should be number'
        gp_parameters['resampling.iterations'] = str(user_input['cc_resampling_iter'])
    if 'cc_seed_value' in user_input:
        assert module_utils.is_number(
            user_input['cc_seed_value']), 'the cc_seed_value should be number'
        gp_parameters['seed.value'] = str(user_input['cc_seed_value'])
    if 'cc_decent_iter' in user_input:
        assert module_utils.is_number(
            user_input['cc_decent_iter']), 'the cc_decent_iter should be number'
        gp_parameters['descent.iterations'] = str(user_input['cc_decent_iter'])
    if 'cc_norm_iter' in user_input:
        assert module_utils.is_number(
            user_input['cc_norm_iter']), 'the cc_norm_iter should be number'
        gp_parameters['normalization.iterations'] = str(user_input['cc_norm_iter'])
    if 'cc_heatmap_size' in user_input:
        assert module_utils.is_number(
            user_input['cc_heatmap_size']), 'the cc_heatmap_size should be number'
        gp_parameters['heat.map.size'] = str(user_input['cc_heatmap_size'])
        
    
        
    gp_parameters['clustering.algorithm'] = parameters['Consensus_algorithm'].upper()
    gp_parameters['distance.measure'] = parameters['cc_distance'].upper()
    gp_parameters['merge.type'] = parameters['merge_type']
    gp_parameters['create.heat.map'] = parameters['create_heatmap']
    gp_parameters['normalize.type'] = parameters['normalize_type']
    gp_parameters['cluster.by'] = parameters['clusterby']
    gp_parameters['resample'] = parameters['cc_resample']

    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module, 'cannot find the %s' % gp_path
    download_directory = outfile
    command = [gp_module, module_name, '-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters', key + ':' + gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory), (
        'there is no output directory for consensusClustering')
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files, 'gene_pattern get error'
    assert module_utils.exists_nz(outfile), (
        'the output file %s for consensusClustering fails' % outfile)
    out_node = bie3.Data(rulebase.ConsensusClusteringFolder,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'ConsensusClusteringFolder' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)
