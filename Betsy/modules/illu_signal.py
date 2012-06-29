#illu_signal.py

import module_utils
import shutil
import os
import rule_engine

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    result_files = os.listdir(single_object.identifier)
    for result_file in result_files:
        if '-controls' not in result_file:
            goal_file = os.path.join(single_object.identifier,result_file)
            shutil.copyfile(goal_file,outfile)
    assert module_utils.exists_nz(outfile),'the output file %s\
                         for illu_signal fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                    parameters,single_object,pipeline,outfile)
    return new_objects



def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_illumina_'+original_file+'.gct'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'illu_folder','contents')
    assert os.path.exists(single_object.identifier),'the input file %s \
             for illu_signal does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
