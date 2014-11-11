#calc_diffexp_with_ttest.py

import os
from Betsy import  module_utils, bie3,rulebase
from genomicode import config
import subprocess


def run(in_nodes, parameters, user_input, network,num_cores):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    diffexp_bin = config.find_diffexp_genes
    assert os.path.exists(diffexp_bin)
    cmd = ['python', diffexp_bin, data_node.identifier,'--cls_file',
           cls_node.identifier, '--algorithm','ttest']
    if 'diffexp_foldchange_value' in user_input:
        foldchange = float(user_input['diffexp_foldchange_value'])
        cmd = cmd + ['--fold_change', str(foldchange)]
    handle = open(outfile,'w')
    try:
        process = subprocess.Popen(cmd, shell=False,
                               stdout=handle,
                               stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    finally:
        handle.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for calc_diffexp_with_ttest fails' % outfile)
    out_node = bie3.Data(rulebase.DiffExprFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 't_test_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


