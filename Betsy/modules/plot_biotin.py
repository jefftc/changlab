#plot_biotin.py
import os
import module_utils
import shutil

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    module_utils.plot_R(identifier,['biotin','housekeeping'],outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

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
        outfile,'biotin_plot',parameters,objects,single_object)
    return new_objects
