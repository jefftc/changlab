#classify_with_svm.py
import svmutil
import sys
import arrayio
import os
from Betsy import read_label_file, module_utils
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    svm_model, data_node_test, cls_node_train = antecedents
    outfile = name_outfile(antecedents, user_options)
    a, train_label, second_line = read_label_file.read(
        cls_node_train.identifier)
    M = arrayio.read(data_node_test.identifier)
    # convert to the format libsvm accept
    test = M.matrix(None, range(len(train_label), M.dim()[1]))
    x_test = module_utils.format_convert(test)
    model = svmutil.svm_load_model(svm_model.identifier)
    a, train_label, second_line = read_label_file.read(
        cls_node_train.identifier)
    y_test = [0] * len(x_test)
    p_label, p_acc, p_val = svmutil.svm_predict(y_test, x_test, model)
    prediction_index = [int(i) for i in p_label]
    prediction = [second_line[i] for i in prediction_index]
    name = test._col_names.keys()[0]
    sample_name = test._col_names[name]
    result = [['Sample_name', 'Predicted_class', 'Confidence']]
    for i in range(len(sample_name)):
        result.append([str(sample_name[i]), prediction[i], str(p_val[i][0])])
    f = file(outfile, 'w')
    for i in result:
        f.write('\t'.join(i))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for classify_with_svm fails' % outfile
    )
    out_node = bie3.Data(rulebase.ClassifyFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    svm_model, data_node_test, cls_node_train = antecedents
    original_file = module_utils.get_inputid(svm_model.identifier)
    filename = 'svm_result' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    svm_model, data_node_test, cls_node_train = antecedents
    identifier = svm_model.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    svm_model_node = module_utils.get_identifier(network, module_id,
                                                 pool, user_attributes,
                                                 datatype='SvmModel')
    data_node_test = module_utils.get_identifier(network, module_id,
                                                 pool, user_attributes,
                                                 contents='class0,class1,test',
                                                 datatype='SignalFile')
    cls_node_train = module_utils.get_identifier(network, module_id,
                                                 pool, user_attributes,
                                                 contents='class0,class1',
                                                 datatype='ClassLabelFile')
    return svm_model_node, data_node_test, cls_node_train
