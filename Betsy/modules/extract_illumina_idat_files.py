#extract_illumina_idat_files.py

import module_utils
import shutil
import os
import Betsy_config
import zipfile
import hash_method
import rule_engine

def extract_all(zipName):
    z = zipfile.ZipFile(zipName)
    for f in z.namelist():
        if f.endswith('/'):
            os.makedirs(f)
        else:
            z.extract(f)

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if zipfile.is_zipfile(identifier):
        directory = os.path.split(identifier)[-1]
        directory = os.path.splitext(directory)[0]
        directory = os.path.join(os.getcwd(),directory)
        extract_all(identifier)
    else:
        directory = identifier
    illumina_file = []
    for filename in os.listdir(directory):
        if filename in ['.DS_Store','._.DS_Store','.Rapp.history']:
            continue
        if filename.endswith('.idat'):
            illumina_file.append(filename)
    if illumina_file:
        os.mkdir(outfile)
        for filename in illumina_file:
            old_file = os.path.join(directory,filename)
            new_file = os.path.join(outfile,filename)
            shutil.copyfile(old_file,new_file)
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
        return new_objects
    else:
        return None
    
def make_unique_hash(identifier,pipeline,parameters):
    original_file = module_utils.get_inputid(identifier)
    hash_profile={'version': 'illumina',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(
        original_file,pipeline,**hash_profile)
    return hash_result

def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'geo_dataset',parameters,objects,single_object)
    return new_objects
    

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

