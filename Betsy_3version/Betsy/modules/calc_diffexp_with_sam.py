#calc_diffexp_with_sam.py
import subprocess
import shutil
import os
import arrayio
from Betsy import module_utils, read_label_file
from Betsy import bie3, rulebase


def run(in_nodes, parameters, user_input, network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    label, label_line, second_line = read_label_file.read(
        cls_node.identifier)
    class_num = len(label)
    assert class_num == 2, (
        'the number of class in %s is not 2' % cls_node.identifier)
    delta = float(parameters['sam_delta'])
    foldchange = float(parameters['sam_foldchange'])
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    sam_script = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'sam_script.py')
    cmd = ['python', sam_script, data_node.identifier,
           cls_node.identifier, outfile, str(delta), str(foldchange)]
    process = subprocess.Popen(cmd, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    process.wait()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for calc_diffexp_with_sam fails' % outfile)
    out_node = bie3.Data(rulebase.DiffExprFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'sam_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


