#clustering.py

import os
import subprocess
import module_utils
import rule_engine
def run(parameters,objects,pipeline):
    """clustering the input file"""
    CLUSTER_BIN = 'cluster'
    distance_para = {'correlation':'1','euclidean':'7'}
    if not int(parameters['k']):
        raise ValueError('K must be an integer')
    dist = distance_para[parameters['distance']]
    alg = {'kmeans': ["-g",dist,"-k",parameters['k']],'pca':["-g",dist,'-pg'],
           'hierarchical':['-m','m'],'som':["-g",dist,'-s']}
    try:    
        com_parameter = alg[parameters['cluster_alg']]
    except:
        raise ValueError("cluster algorithm is not recognized")
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    command = [CLUSTER_BIN,'-f',identifier,'-u',outfile]
    for i in com_parameter:
        command.append(i)
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    result_files = os.listdir(os.getcwd())
    fullpath = [os.path.realpath(i) for i in result_files]
    sizelist = [os.path.getsize(i) for i in fullpath]
    min_size = min(sizelist)
    min_index = sizelist.index(min_size)
    outputfile = result_files[min_index]
    os.rename(outputfile,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')
    
def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','cluster_file')

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
