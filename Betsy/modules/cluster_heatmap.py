#cluster_heatmap.py

import os
import hash_method
import subprocess
import module_utils
import rule_engine
def run(parameters,objects,pipeline):
    """generate a heatmap of input file"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    import Betsy_config
    Heatmap_BIN = Betsy_config.ARRAYPLOT
    color='--color=' + parameters['color'].replace('_','-')
    command = ['python', Heatmap_BIN,
                   "-y","1","-x","200","--label_arrays","--no_autoscale",
               color,identifier,'-o',outfile]
    
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
            identifier,pipeline,parameters)
   
def get_identifier(parameters,objects):
    if parameters['cluster_alg'] == 'no_cluster_alg':
        identifier,single_object = module_utils.find_object(
            parameters,objects,'signal_file','Contents,DatasetId')
    else:
        identifier,single_object = module_utils.find_object(
            parameters,objects,'cluster_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    cwd = os.getcwd()
    outfile = os.path.join(cwd,os.path.split(identifier)[-1]+'.png')
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'cluster_heatmap',parameters,objects,single_object)
    return new_objects
    
