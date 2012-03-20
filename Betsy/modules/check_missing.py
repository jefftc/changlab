#check_missing.py
import os
import shutil
import module_utils

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if is_missing(parameters,objects):
        shutil.copyfile(identifier,outfile)
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline)
        return new_objects
    else:
        return None
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
    
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
