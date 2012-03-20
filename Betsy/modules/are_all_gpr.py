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
    outfile = get_outfile(parameters,objects,pipeline)
    check = []
    for i in os.listdir(identifier):
        if i == '.DS_Store':
            pass
        else:
            gpr = gpr_module.check_gpr(os.path.join(identifier,i))
            check.append(gpr)
    if check == [True]*len(check):
        os.mkdir(outfile)
        for i in os.listdir(identifier):
            if i != '.DS_Store':
                old_file = os.path.join(identifier,i)
                new_file = os.path.join(outfile,i)
                shutil.copyfile(old_file,new_file)
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline)
        return new_objects
    else:
         return None


def make_unique_hash(identifier,pipeline,parameters):
    inputid = module_utils.get_inputid(identifier)
    hash_profile = {'version': 'gpr',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result = hash_method.hash_parameters(
                 inputid,pipeline,**hash_profile)
    return hash_result


def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_'+ hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        ooutfile,'geo_dataset',parameters,objects,single_object)
    return new_objects


def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'geo_dataset','Contents,DatasetId')
    assert os.path.exists(identifier),'folder %s does not exit.' % identifier
    assert os.path.isdir(identifier),"input is not a folder"
    return identifier,single_object

 
