#annotate_genes_with_gsea.py
import subprocess
import shutil
import os
import arrayio
from genomicode import arrayplatformlib, config
from Betsy import module_utils, read_label_file
from Betsy import bie3
from Betsy import rulebase

def run(in_nodes, parameters,user_input, network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    gsea_path = config.gsea
    gsea_module = module_utils.which(gsea_path)
    assert gsea_module, 'cannot find the %s' % gsea_path
    M = arrayio.read(data_node.identifier)
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    chipname = x[0][1]
    platform = chipname + '.chip'
    download_directory = os.path.join(os.getcwd(),'gsea_result')
    command = [gsea_module, data_node.identifier, '--cls_file',
               cls_node.identifier, '--platform', platform,
               download_directory]
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(download_directory), (
        'there is no output directory for GSEA')
    shutil.copytree(download_directory, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annotate_genes_with_gsea fails' % outfile)
    out_node = bie3.Data(rulebase.GseaFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='PrettySignalFile')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'gsea_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)



