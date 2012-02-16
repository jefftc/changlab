#cluster_heatmap.py

import os
import hash_method
import subprocess
import module_utils
import rule_engine
def run(parameters,objects):
    """generate a heatmap of input file"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects =  get_outfile(parameters,objects)
    import Betsy_config
    Heatmap_BIN = Betsy_config.ARRAYPLOT
    color='--color='+parameters['color'].replace('_','-')
    command = ['python', Heatmap_BIN,
                   "-y","1","-x","200","--label_arrays",
               "--no_autoscale",color,identifier]
    
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    outputfile = os.listdir(os.getcwd())[0]
    os.rename(outputfile,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects
    
def make_unique_hash(parameters,objects):
    if parameters['cluster_alg']=='no_cluster_alg':
        return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')
    else:
        return module_utils.make_unique_hash(parameters,objects,'cluster_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    if parameters['cluster_alg']=='no_cluster_alg':
        return module_utils.get_outfile(parameters,objects,
                                    'signal_file','Contents,DatasetId','cluster_heatmap')
    else:
        return module_utils.get_outfile(parameters,objects,
                                    'cluster_file','Contents,DatasetId','cluster_heatmap')
def get_identifier(parameters,objects):
    if parameters['cluster_alg']=='no_cluster_alg':
        return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
    else:
        return module_utils.find_object(parameters,objects,'cluster_file','Contents,DatasetId')
