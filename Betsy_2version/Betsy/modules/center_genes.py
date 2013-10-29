#center_genes.py
import os
import subprocess
from Betsy import module_utils,bie,rulebase


def run(data_node,parameters, network):
    """mean or median"""
    CLUSTER_BIN = 'cluster'
    center_alg = {'mean': 'a', 'median': 'm'}
    try:
        center_parameter = center_alg[parameters['gene_center']]
    except:
        raise ValueError("Centering parameter is not recognized")
    outfile = name_outfile(data_node)
    process = subprocess.Popen([CLUSTER_BIN, '-f', data_node.attributes['filename'],
                                '-cg', center_parameter, '-u', outfile],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    outputfile = outfile + '.nrm'
    os.rename(outputfile, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for centering fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile1,**new_parameters)
    return out_node


def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_center_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    new_parameters['format']='tdf'
    return new_parameters
    

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
