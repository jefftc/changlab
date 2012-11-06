#apply_svm_model.py
import svmutil
import sys
import arrayio
import os
from Betsy import read_label_file, module_utils


def run(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    model_file = module_utils.find_object(
        parameters, objects, 'svm_model', 'traincontents')
    assert os.path.exists(model_file.identifier), (
        'the model file %s for svm test does not exist'
        % model_file.identifier)
    test = arrayio.read(single_object.identifier)
    # convert to the format libsvm accept
    x_test = module_utils.format_convert(test)
    test_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'testcontents')
    model = svmutil.svm_load_model(model_file.identifier)
    if test_label_file:
        assert os.path.exists(test_label_file.identifier), (
            'the test_label_file %s for svm test does not exist'
            % test_label_file.identifier)
        a, test_label, second_line = read_label_file.read(
            test_label_file.identifier)
        actual_label = [second_line[int(i)] for i in test_label]
        y_test = [int(x) for x in test_label]
    else:
        y_test = [0] * test.ncol()
    p_label, p_acc, p_val = svmutil.svm_predict(y_test, x_test, model)
    train_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'traincontents')
    assert os.path.exists(train_label_file.identifier), (
        'the train_label_file %s for svm test does not exist'
        % train_label_file.identifier)
    result, label_line, second_line = read_label_file.read(
        train_label_file.identifier)
    prediction_index = [int(i) for i in p_label]
    prediction = [second_line[i] for i in prediction_index]
    name = test._col_names.keys()[0]
    sample_name = test._col_names[name]
    if test_label_file:
        result = [['Sample_name', 'Predicted_class', 'Confidence',
                   'Actual_class', 'Correct?']]
        for i in range(len(sample_name)):
            if prediction[i] == actual_label[i]:
                correct = 'yes'
            else:
                correct = 'no'
            result.append(
                [sample_name[i], prediction[i], str(p_val[i][0]),
                 actual_label[i], correct])
    else:
        result = [['Sample_name', 'Predicted_class', 'Confidence']]
        for i in range(len(sample_name)):
            result.append(
                [sample_name[i], prediction[i], str(p_val[i][0])])
    f = file(outfile, 'w')
    for i in result:
        f.write('\t'.join(i))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for svm test fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'prediction_result_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'testcontents')
    assert os.path.exists(single_object.identifier), (
        'the test file %s for svm test does not exist'
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
