#check_illumina.py

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

def run(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    if zipfile.is_zipfile(identifier):
        directory = os.path.split(identifier)[-1]
        directory = os.path.splitext(directory)[0]
        directory = os.path.join(os.getcwd(),directory)
        extract_all(identifier)
    else:
        directory = identifier
    flag = True
    for filename in os.listdir(directory):
        if filename in ['.DS_Store','._.DS_Store']:
            pass
        else:
            if not filename.endswith('.idat'):
                flag = False
    if flag:
        os.mkdir(outfile)
        for filename in os.listdir(directory):
            if filename not in ['.DS_Store','._.DS_Store']:
                old_file = os.path.join(directory,filename)
                new_file = os.path.join(outfile,filename)
                shutil.copyfile(old_file,new_file)
        
        module_utils.write_Betsy_parameters_file(parameters,single_object)
        return new_objects
    else:
        return None
    
def make_unique_hash(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    hash_profile={'filename':os.path.split(identifier)[-1],
                   'version': 'illumina',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(**hash_profile)
    return hash_result

def get_outfile(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    old_filename = os.path.split(identifier)[-1]
    if '_BETSYHASH_' in old_filename: 
        original_file = '_'.join(old_filename.split('_')[:-2])
    else:
        original_file = old_filename
    hash_string = make_unique_hash(parameters,objects)
    filename = original_file +'_BETSYHASH_'+ hash_string
    outfile = os.path.join(os.getcwd(),filename)
    attributes = parameters.values()
    objecttype = 'geo_dataset'
    new_object = rule_engine.DataObject(objecttype,attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
    
