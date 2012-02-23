#are_all_gpr.py
import os
import gzip
import shutil
import rule_engine
import module_utils
import hash_method
import gpr_module

def run(parameters,objects,pipeline):
    """check if all the file are gpr format"""
    identifier,single_object = get_identifier(parameters,objects)
    assert os.path.exists(identifier),'folder %s does not exit.' % identifier
    assert os.path.isdir(identifier),"input is not a folder"
    outfile,new_objects = get_outfile(parameters,objects)
    check=[]
    for i in os.listdir(identifier):
        if i == '.DS_Store':
            pass
        else:
            gpr = gpr_module.check_gpr(os.path.join(identifier,i))
            check.append(gpr)
    if check==[True]*len(check):
        os.mkdir(outfile)
        for i in os.listdir(identifier):
            if i != '.DS_Store':
                old_file = os.path.join(identifier,i)
                new_file = os.path.join(outfile,i)
                shutil.copyfile(old_file,new_file)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
        return new_objects
    else:
         return None


def make_unique_hash(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    hash_profile={'filename':os.path.split(identifier)[-1],
                   'version': 'gpr',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(**hash_profile)
    return hash_result


def get_outfile(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    old_filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in old_filename: 
        original_file = '_'.join(old_filename.split('_')[:-2])
    else:
        original_file = old_filename
    hash_string = make_unique_hash(parameters,objects)
    filename = original_file + hash_string
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


