#relabel_samples.py
import os
import shutil
import arrayio
import subprocess
from genomicode import config
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(in_nodes,parameters, network):
    data_node,rename_node = in_nodes
    outfile = name_outfile(in_nodes)
    rename_path = config.slice_matrix
    rename_BIN = module_utils.which(rename_path)
    assert rename_BIN,'cannot find the %s' %rename_path
    command=['python',rename_BIN,data_node.attributes['filename'],
             '--relabel_col_ids',
             rename_node.attributes['filename']+',NewName']
    f=file(outfile,'w')
    process = subprocess.Popen(command,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
    f.close()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
     
    assert module_utils.exists_nz(outfile),(
        'the output file %s for relabel_samples does not exist'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node



def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile')
    rename_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='RenameFile')
    return data_node, rename_node

def name_outfile(in_nodes):
    data_node, rename_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_rename_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,rename_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

