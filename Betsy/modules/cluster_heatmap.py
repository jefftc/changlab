#cluster_heatmap.py

import os
import hash_method
import subprocess
import module_utils
import rule_engine
def run(parameters,objects,pipeline,options=None):
    """generate a heatmap of input file"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    import Betsy_config
    Heatmap_path = Betsy_config.ARRAYPLOT
    Heatmap_BIN = module_utils.which(Heatmap_path)
    assert Heatmap_BIN,'cannot find the %s' %Heatmap_path
    command = ['python', Heatmap_BIN,identifier,'-o',outfile,"--label_arrays",
               "--grid"]
    if 'color' in parameters.keys():
        color=['--color' , parameters['color'].replace('_','-')]
        command.extend(color)
    if 'hm_width' in parameters.keys():
        command.extend(['-x',parameters['hm_width']])
    if 'hm_height' in parameters.keys():
        command.extend(['-y',parameters['hm_height']])
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),'the output file %s\
             for cluster_heatmap fails' %outfile
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
            parameters,objects,'signal_file','contents')
    else:
        identifier,single_object = module_utils.find_object(
            parameters,objects,'cluster_file','contents')
    assert os.path.exists(identifier),'the input file %s\
            for cluster_heatmap does not exist'%identifier
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
    
