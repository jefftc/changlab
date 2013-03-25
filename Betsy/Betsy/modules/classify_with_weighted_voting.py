#classify_with_weighted_voting.py
import shutil
import os
import subprocess
from genomicode import config
from Betsy import module_utils, read_label_file
from time import strftime,localtime

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    train_identifier = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    train_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'traincontents')
    test_label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'testcontents')
    test_file = module_utils.find_object(
        parameters, objects, 'signal_file', 'testcontents')
    if 'given' in test_label_file.attributes:
        actual = True
    else:
        actual = False
    assert os.path.exists(test_file.identifier), (
        'the test file %s for weighted_voting does not exist'
        % test_file.identifier)
    assert os.path.exists(test_label_file.identifier), (
        'cannot find test_label_file %s for weighted_voting'
        % test_label_file.identifier)
    assert os.path.exists(train_label_file.identifier), (
        'cannot find train_label_file %s for weighted_voting'
        % train_label_file.identifier)
    module_name = 'WeightedVoting'
    gp_parameters = dict()
    file1, file2 = module_utils.convert_to_same_platform(
        train_identifier.identifier, test_file.identifier)
    gp_parameters['train.filename'] = file1
    gp_parameters['train.class.filename'] = train_label_file.identifier
    gp_parameters['test.filename'] = file2
    gp_parameters['test.class.filename'] = test_label_file.identifier
    num_features = parameters['num_features']
    if int(num_features) > 0:
        gp_parameters['num.features'] = str(parameters['num_features'])
    if 'wv_minstd' in parameters.keys():
        assert module_utils.is_number(
            parameters['wv_minstd']), 'the sv_minstd should be number'
        gp_parameters['min.std'] = str(parameters['wv_minstd'])
    wv_feature_stat = ['wv_snr', 'wv_ttest', 'wv_snr_median',
                       'wv_ttest_median', 'wv_snr_minstd',
                       'wv_ttest_minstd', 'wv_snr_median_minstd',
                       'wv_ttest_median_minstd']
    if 'wv_feature_stat' in parameters.keys():
        assert parameters['wv_feature_stat'] in wv_feture_stat, (
            'the wv_feature_stat is invalid')
        gp_parameters['feature.selection.statistic'] = str(
            wv_feature_stat.index(parameters[
                'wv_feature_stat']))
    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module, 'cannot find the %s' % gp_path
    download_directory = os.path.join(os.getcwd(),'wv_result')
    command = [gp_module, module_name, '-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters', key + ':' + gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory), (
        'there is no output directory for weightedVoting')
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files, 'gene_pattern get error'
    gp_files = os.listdir(download_directory)
    result, label_line, second_line = read_label_file.read(
        train_label_file.identifier)
    for gp_file in gp_files:
        if gp_file.endswith('pred.odf'):
            gp_file = os.path.join(download_directory, gp_file)
            f = file(gp_file, 'r')
            text = f.readlines()
            f.close()
            os.rename(os.path.join(download_directory, gp_file),
                      os.path.join(download_directory, 'prediction.odf'))
            assert text[1][0: 12] == 'HeaderLines='
            start = int(text[1][12: -1])
            if actual:
                newresult = [['Sample_name',
                              'Predicted_class', 'Confidence',
                              'Actual_class', 'Correct?']]
            else:
                newresult = [['Sample_name', 'Predicted_class', 'Confidence']]
            for i in text[start + 2:]:
                line = i.split()
                n = len(line)
                if actual:
                    if line[n - 3] == line[n - 4]:
                        correct = 'yes'
                    else:
                        correct = 'no'
                    newline = [' '.join(line[0: n - 4]), line[n - 3],
                               line[n - 2], line[n - 4], correct]
                else:
                    newline = [' '.join(line[0: n - 4]), line[n - 3],
                               line[n - 2]]
                newresult.append(newline)
            f = file(outfile, 'w')
            for i in newresult:
                f.write('\t'.join(i))
                f.write('\n')
            f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for classify_with_weighted_voting fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, [train_identifier,test_file,train_label_file,test_label_file],
        pipeline, outfile,starttime,user,jobname)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'traincontents')
    assert os.path.exists(single_object.identifier), (
        'the train file %s for classify_with_weighted_voting does not exist'
        % single_object.identifier)
    return single_object


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'weighted_voting_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


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
