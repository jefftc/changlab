#classify_with_random_forest.py

import sys
import arrayio
import os
from Betsy import read_label_file, module_utils
from genomicode import jmath


def run(parameters, objects, pipeline):
    train_file = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    train_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'traincontents')
    test_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'testcontents')
    test_file = module_utils.find_object(
        parameters, objects, 'signal_file', 'testcontents')
    assert os.path.exists(train_label_file.identifier), (
        'cannot find train_label_file %s for classify_with_random_forest'
        % train_label_file.identifier)
    assert os.path.exists(test_file.identifier), (
        'the test file %s for classify_with_random_forest does not exist'
        % test_file.identifier)
    result, label_line, second_line = read_label_file.read(
        train_label_file.identifier)
    y = ['"' + second_line[int(i)] + '"' for i in label_line]
    if test_label_file:
        assert os.path.exists(test_label_file.identifier), (
            'test_label_file %s for classify_with_random_forest does not exist'
            % test_label_file.identifier)
        a, test_label, second_line = read_label_file.read(
            test_label_file.identifier)
        actual_label = [second_line[int(i)] for i in test_label]
    R = jmath.start_R()
    M_train = arrayio.read(train_file.identifier)
    M1 = M_train.slice()
    M_train = jmath.transpose(M1)
    jmath.R_equals_matrix(M_train, 'data')
    M_test = arrayio.read(test_file.identifier)
    M2 = M_test.slice()
    M2 = jmath.transpose(M2)
    jmath.R_equals_matrix(M2, 'test')
    jmath.R_equals(y, 'y')
    R('y<-as.factor(y)')
    R('require(randomForest,quietly=TRUE)')
    R('library(randomForest)')
    R('model<-randomForest(data,y=y,importance=TRUE)')
    R('predict_result<-predict(model, test)')
    predict_result = R['predict_result']
    levels = predict_result.levels
    predict_labels = predict_result[:]
    predict_labels = [levels[i - 1] for i in predict_labels]
    name = M_test._col_names.keys()[0]
    sample_name = M_test._col_names[name]
    if test_label_file:
        result = [['Sample_name', 'Predicted_class',
                   'Confidence', 'Actual_class', 'Correct?']]
        for i in range(len(sample_name)):
            if predict_labels[i] == actual_label[i]:
                correct = 'yes'
            else:
                correct = 'no'
            result.append(
                [sample_name[i], predict_labels[i], '',
                 actual_label[i], correct])
    else:
        result = [['Sample_name', 'Predicted_class', 'Confidence']]
        for i in range(len(sample_name)):
            result.append(
                [sample_name[i], predict_labels[i], ''])
    f = file(outfile, 'w')
    for i in result:
        f.write('\t'.join(i))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for classify_with_random_forest fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, train_file, pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'random_forest_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'traincontents')
    assert os.path.exists(single_object.identifier), (
        'the test file %s for classify_with_random_forest does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    parameters = module_utils.renew_parameters(parameters, ['status'])
    attributes = parameters.values()
    new_object = module_utils.DataObject(
        'classification_file', attributes, outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
