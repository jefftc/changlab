#normalize_samples_with_bfrm.py
import os
import subprocess
import shutil
import arrayio
import tempfile
from genomicode import config
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, network):
    outfile = name_outfile(data_node)
    bfrm_path = config.bfrmnorm
    bfrm_BIN = module_utils.which(bfrm_path)
    assert bfrm_BIN,'cannot find the %s' %bfrm_path
    num_factor = 1
    #num_factor = 10
    if 'num_factors' in parameters.keys():
        num_factor = int(parameters['num_factors'])
        assert num_factor >= 1, 'the num_factor should be >=1'
        M = arrayio.read(single_object.identifier)
        col_num = M.ncol()
        assert num_factor <= col_num,(
            'the num_factor should be less than %d'%col_num)
    tmp = 'tmp_dir'
    command = ['python', bfrm_BIN,data_node.attributes['filename'],'-f',str(num_factor), '-o',tmp]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tmp),(
        'the output dir %s for bfrm_normalize fails' %tmp)
    assert module_utils.exists_nz(os.path.join(tmp,'normalized.gct')),(
        'the output gct file for bfrm_normalize fails')
    out = os.path.join(tmp,'normalized.gct')
    M = arrayio.read(out)
    M_new = arrayio.convert(M,to_format = arrayio.pcl_format)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(SignalFile_rule.SignalFile, **new_parameters)
    return out_node

    

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_bfrm_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)