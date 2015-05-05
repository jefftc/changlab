#run_loocv.py
import arrayio
import os
import svmutil
from Betsy import bie3, rule_engine_bie3
from Betsy import rulebase
from Betsy import read_label_file
from Betsy import module_utils
from genomicode import config
import subprocess


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node_train, cls_node_train = antecedents
    outfile = name_outfile(antecedents, user_options)
    module_name = 'WeightedVotingXValidation'
    module_id_version = '00028:2'
    gp_parameters = dict()
    file1 = data_node_train.identifier
    result, label_line, class_name = read_label_file.read(
        cls_node_train.identifier)
    gp_parameters['data.filename'] = file1
    gp_parameters['class.filename'] = cls_node_train.identifier
    if 'wv_num_features' in user_options:
        gp_parameters['num.features'] = str(user_options['wv_num_features'])
##    if 'wv_minstd' in user_input:	
##    	assert module_utils.is_number(
##            user_input['wv_minstd']), 'the sv_minstd should be number'
##        gp_parameters['min.std'] = str(user_input['wv_minstd'])
##        
##    wv_feature_stat = ['wv_snr', 'wv_ttest', 'wv_snr_median',
##                       'wv_ttest_median', 'wv_snr_minstd',
##                       'wv_ttest_minstd', 'wv_snr_median_minstd',
##                       'wv_ttest_median_minstd']
##    
##    assert parameters['wv_feature_stat'] in wv_feature_stat, (
##            'the wv_feature_stat is invalid')
##    gp_parameters['feature.selection.statistic'] = str(
##            wv_feature_stat.index(parameters[
##                'wv_feature_stat']))

    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module, 'cannot find the %s' % gp_path
    download_directory = os.path.join(os.getcwd(), 'wv_result')
    command = [gp_module, module_name, '--id_and_version', module_id_version,
               '-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters', key + ':' + gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory), (
        'there is no output directory for weightedVotingXValidation'
    )
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files, 'gene_pattern get error'
    gp_files = os.listdir(download_directory)
    for gp_file in gp_files:
        if gp_file.endswith('pred.odf'):
            gp_file = os.path.join(download_directory, gp_file)
            f = file(gp_file, 'r')
            text = f.readlines()
            f.close()
            os.rename(os.path.join(download_directory, gp_file),
                      os.path.join(download_directory, 'prediction.odf'))
            assert text[1][0:12] == 'HeaderLines='
            start = int(text[1][12:-1])
            newresult = [['Sample_name', 'Predicted_class', 'Confidence',
                          'Actual_class', 'Correct?']]
            for i in text[start + 2:]:
                line = i.split()
                n = len(line)
                newline = [' '.join(line[0:n - 4]), line[n - 3], line[n - 2],
                           line[n - 4], line[n - 1]]
                newresult.append(newline)
            f = file(outfile, 'w')
            for i in newresult:
                f.write('\t'.join(i))
                f.write('\n')
            f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for run_loocv_weighted_voting fails' % outfile
    )
    out_node = bie3.Data(rulebase.ClassifyFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(
        datatype_name='SignalFile', contents="class0,class1")
    filter2 = module_utils.AntecedentFilter(
        datatype_name='ClassLabelFile', contents="class0,class1")
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'predication_loocv_wv' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
