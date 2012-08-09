#weighted_voting.py
import module_utils
import shutil
import os
from genomicode import jmath
import Betsy_config
import subprocess
import read_label_file
def run(parameters,objects,pipeline):
    train_identifier = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    train_label_file= module_utils.find_object(
        parameters,objects,'class_label_file','traincontents')
    test_label_file = module_utils.find_object(
        parameters,objects,'class_label_file','testcontents')
    test_file = module_utils.find_object(
        parameters,objects,'signal_file','testcontents')
    if 'given' in test_label_file.attributes:
        actual = True
    else:
        actual = False
    assert os.path.exists(test_file.identifier),(
        'the test file %s for weighted_voting does not exist'
        %test_file.identifier)
    assert os.path.exists(test_label_file.identifier),(
        'cannot find test_label_file %s for weighted_voting'
        %test_label_file.identifier)
    assert os.path.exists(train_label_file.identifier),(
        'cannot find train_label_file %s for weighted_voting'
        %train_label_file.identifier)
    module_name = 'WeightedVoting'
    gp_parameters = dict()
    gp_parameters['train.filename'] = train_identifier.identifier
    gp_parameters['train.class.filename'] = train_label_file.identifier
    gp_parameters['test.filename'] = test_file.identifier
    gp_parameters['test.class.filename'] = test_label_file.identifier
    if 'wv_num_features' in parameters.keys():
        assert parameters['wv_num_features'].isdigit(
            ),'the wv_num_features should be digit numbers'
        gp_parameters['num.features'] = str(parameters['wv_num_features'])
    if 'wv_minstd' in parameters.keys():
        assert module_utils.is_number(
            parameters['wv_minstd']), 'the sv_minstd should be number'
        gp_parameters['min.std'] = str(parameters['wv_minstd'])
    wv_feature_stat = ['wv_snr','wv_ttest','wv_snr_median','wv_ttest_median',
                       'wv_snr_minstd','wv_ttest_minstd','wv_snr_median_minstd',
                       'wv_ttest_median_minstd']
    if 'wv_feature_stat' in parameters.keys():
        assert parameters['wv_feature_stat'] in wv_feture_stat, 'the wv_feature_stat is invalid'
        gp_parameters['feature.selection.statistic'] = str(
                               wv_feature_stat.index(parameters['wv_feature_stat']))
    gp_path = Betsy_config.GENEPATTERN
    gp_module = module_utils.which(gp_path)
    assert gp_module,'cannot find the %s' %gp_path
    command = [gp_module, module_name]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    out_text =  process.stdout.read()
    out_lines = out_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory),(
        'there is no output directory for weightedVoting')
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    gp_files = os.listdir(download_directory)
    result,label_line,second_line = read_label_file.read(train_label_file.identifier)
    for gp_file in gp_files:
        if gp_file.endswith('pred.odf'):
            gp_file = os.path.join(download_directory,gp_file)
            f = file(gp_file,'r')
            text = f.readlines()
            assert text[1][0:12]=='HeaderLines='
            start=int(text[1][12:-1])
            newresult=[['Sample_name','Predicted_class','Confidence','Actual_class']]
            for i in text[start+2:]:
                line = i.split()
                n = len(line)
                if actual:
                    newline = [' '.join(line[0:n-4]),line[n-3],line[n-2],line[n-4]]
                else:
                    newline = [' '.join(line[0:n-4]),line[n-3],line[n-2],'']
                newresult.append(newline)
            f=file(outfile,'w')
            for i in newresult:
                f.write('\t'.join(i))
                f.write('\n')
            f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for weighted_voting fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,train_identifier,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','traincontents')
    assert os.path.exists(single_object.identifier),(
        'the train file %s for weighted_voting does not exist'
        %single_object.identifier)
    return single_object

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'weighted_voting_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'weighted_voting',parameters,objects,single_object)
    return new_objects
