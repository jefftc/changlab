#weighted_voting.py
import module_utils
import shutil
import os
from genomicode import jmath
import Betsy_config

def run(parameters,objects,pipeline):
    train_identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    train_label_file,obj = module_utils.find_object(parameters,objects,'class_label_file','TrainContents,Train_DatasetId')
    test_label_file,obj = module_utils.find_object(parameters,objects,'class_label_file','TestContents,Test_DatasetId')
    test_file = module_utils.find_object(parameters,objects,'signal_file','TestContents,Test_DatasetId')
    assert os.path.exists(label_file),'cannot find label_file'
    module_name = 'WeightedVoting'
    parameters = dict()
    parameters['train.filename'] = train_identifier
    parameters['train.class.filename'] = train_label_file
    parameters['test.filename'] = test_file
    parameters['test.class.filename'] = test_label_file
    gp_module = Betsy_config.GENEPATTERN
    command = [gp_module, module_name]
    for key in parameters.keys():
        a = ['--parameters',key+':'+ parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    output_text =  process.stdout.read()
    out_lines = output_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    os.rename(download_directory,outfile)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','TrainContents,Train_DatasetId',pipeline)

def get_identifier(parameters,objects,pipeline):
    return module_utils.find_object(
        parameters,objects,'signal_file','TrainContents,Train_DatasetId',pipeline)

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','TrainContents,Train_DatasetId','signal_file')
    
