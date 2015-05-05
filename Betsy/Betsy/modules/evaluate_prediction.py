#evaluate_prediciton

import os
from genomicode import config
from Betsy import module_utils, read_label_file
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node_test = antecedents
    outfile = name_outfile(antecedents, user_options)
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
    out_node = bie3.Data(rulebase.ClassifyFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    data_node, cls_node_test = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'prediction_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node_test = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='ClassifyFile')
    filter2 = module_utils.AntecedentFilter(
        datatype_name='ClassLabelFile', contents="test")
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x
