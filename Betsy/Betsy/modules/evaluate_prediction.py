#evaluate_prediciton

import shutil
import os
import subprocess
from genomicode import config
from Betsy import module_utils, read_label_file
from Betsy import bie3
from Betsy import rulebase


def run(in_nodes, parameters, user_input, network, num_cores):
    data_node, cls_node_test = in_nodes
    outfile = name_outfile(in_nodes, user_input)
    f = file(data_node.identifier, 'r')
    text = f.readlines()
    f.close()
    a, test_label, second_line = read_label_file.read(cls_node_test.identifier)
    actual_label = [second_line[int(i)] for i in test_label]
    f = file(outfile, 'w')
    header = text[0].replace('\n', '').split('\t')
    header.extend(['Actual_class', 'Correct?'])
    f.write('\t'.join(header) + '\n')
    for index in range(1, len(text)):
        line = text[index].replace('\n', '')
        line = line.split('\t')
        correct = 'no'
        if line[1] == actual_label[index - 1]:
            correct = 'yes'
        line.extend([actual_label[index - 1], correct])
        f.write('\t'.join(line) + '\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for evaluate_prediction' % outfile
    )
    out_node = bie3.Data(rulebase.ClassifyFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(in_nodes, user_input):
    data_node, cls_node_test = in_nodes
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'prediction_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, in_nodes):
    return parameters


def make_unique_hash(in_nodes, pipeline, parameters, user_input):
    data_node, cls_node_test = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes,
                                            datatype='ClassifyFile')
    cls_node_test = module_utils.get_identifier(network, module_id, data_nodes,
                                                user_attributes,
                                                contents='test',
                                                datatype='ClassLabelFile')
    return data_node, cls_node_test
