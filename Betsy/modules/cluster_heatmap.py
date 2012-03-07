#cluster_heatmap.py

import os
import hash_method
import subprocess
import module_utils
import rule_engine
def run(parameters,objects,pipeline):
    """generate a heatmap of input file"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects =  get_outfile(parameters,objects,pipeline)
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
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(parameters,objects,pipeline):
    if parameters['cluster_alg'] == 'no_cluster_alg':
        return module_utils.make_unique_hash(
            parameters,objects,'signal_file','Contents,DatasetId',pipeline)
    else:
        return module_utils.make_unique_hash(
            parameters,objects,'cluster_file','Contents,DatasetId',pipeline)

def get_identifier(parameters,objects):
    if parameters['cluster_alg'] == 'no_cluster_alg':
        return module_utils.find_object(
            parameters,objects,'signal_file','Contents,DatasetId')
    else:
        return module_utils.find_object(
            parameters,objects,'cluster_file','Contents,DatasetId')
    
def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    cwd = os.getcwd()
    outfile = os.path.join(cwd,os.path.split(identifier)[-1]+'.png')
    if 'Status' in parameters.keys():
        del parameters['Status']
    attributes = parameters.values()
    new_object = rule_engine.DataObject('cluster_heatmap',attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects
