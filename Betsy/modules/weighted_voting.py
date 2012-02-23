#weighted_voting.py
import module_utils
import shutil
import os
from genomicode import jmath
import Betsy_config

def run(parameters,objects,pipeline):
    train_identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    train_label_file,obj=module_utils.find_object(parameters,objects,'class_label_file','TrainContents,Train_DatasetId')
    test_label_file,obj=module_utils.find_object(parameters,objects,'class_label_file','TestContents,Test_DatasetId')
    test_file=module_utils.find_object(parameters,objects,'signal_file','TestContents,Test_DatasetId')
    assert os.path.exists(label_file),'cannot find label_file'
    module_name='WeightedVoting'
    parameters=dict()
    parameters['train.filename']=train_identifier
    parameters['train.class.filename']=train_label_file
    parameters['test.filename']=test_file
    parameters['test.class.filename']=test_label_file
    download_directory=module_utils.run_gp_module(module_name,parameters)
    os.rename(download_directory,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','TrainContents,Train_DatasetId')

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','TrainContents,Train_DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','TrainContents,Train_DatasetId','signal_file')
    
