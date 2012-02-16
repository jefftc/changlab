#class_neighbors.py
import module_utils
import shutil
import os
from genomicode import jmath
import Betsy_config

def run(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    label_file,obj=module_utils.find_object(parameters,objects,'class_label_file','Contents,DatasetId')
    assert os.path.exists(label_file),'cannot find label_file'
    module_name='ClassNeighbors'
    parameters=dict()
    parameters['data.filename']=identifier
    parameters['class.filename']=label_file
    download_directory=module_utils.run_gp_module(module_name,parameters)
    os.rename(download_directory,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
