#annotate_genes_with_gsea.py
import subprocess
import shutil
import os
import arrayio
from genomicode import arrayplatformlib, config
#from Betsy
import module_utils, read_label_file
from time import strftime,localtime
import bie
import rulebase

def run(in_nodes, parameters):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    gsea_path = config.gsea
    gsea_module = module_utils.which(gsea_path)
    assert gsea_module, 'cannot find the %s' % gsea_path
    M = arrayio.read(data_node.attributes['filename'])
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    chipname = x[0][1]
    platform = chipname + '.chip'
    download_directory = os.path.join(os.getcwd(),'gsea_result')
    command = [gsea_module, data_node.attributes['filename'], '--cls_file',
               cls_node.attributes['filename'], '--platform', platform,
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
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.GseaFile,**new_parameters)
    return out_node


def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'gsea_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)



