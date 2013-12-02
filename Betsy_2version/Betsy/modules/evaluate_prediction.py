#evaluate_prediciton

import shutil
import os
import subprocess
from genomicode import config
from Betsy import module_utils, read_label_file
from Betsy import bie 
from Betsy import rulebase

def run(in_nodes,parameters, network):
    data_node,cls_node_test = in_nodes
    outfile = name_outfile(in_nodes)
    f = file(data_node.attributes['filename'],'r')
    text = f.readlines()
    f.close()
    a, test_label, second_line = read_label_file.read(
            cls_node_test.attributes['filename'])
    actual_label = [second_line[int(i)] for i in test_label]
    f =  file(outfile,'w')
    header = text[0].replace('\n','').split('\t')
    header.extend(['Actual_class', 'Correct?'])
    f.write('\t'.join(header)+'\n')
    for index in range(1,len(text)):
        line = text[index].replace('\n','')
        line = line.split('\t')
        correct = 'no'
        if line[1] ==  actual_label[index-1]:
            correct = 'yes'
        line.extend([actual_label[index-1],correct])
        f.write('\t'.join(line)+'\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for evaluate_prediction' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.ClassifyFile,**new_parameters)
    return out_node

def name_outfile(in_nodes):
    data_node,cls_node_test = in_nodes
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'prediction_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node_test = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,
                                            datatype='ClassifyFile')
    cls_node_test = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='test',
                                            datatype='ClassLabelFile')
    return data_node,cls_node_test
