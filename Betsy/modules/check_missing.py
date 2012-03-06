#check_missing.py
import os
import shutil
import module_utils

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    if is_missing(parameters,objects):
        shutil.copyfile(identifier,outfile)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline)
        return new_objects
    else:
        return None
    
def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','signal_file',pipeline)
    
def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')


def is_missing(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    import arrayio
    M = arrayio.read(identifier)
    has_missing = False
    for i in range(M.dim()[0]):
       for j in range(M.dim()[1]):
           if M._X[i][j] == None:
               has_missing = True
               break
       if has_missing:
            break
    return has_missing
