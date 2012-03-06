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
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
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
    result_formats = {'kmeans': 'cdt','pca':'pc.txt',
           'hierarchical':'nrm','som':'txt'}
    result_format = result_formats[parameters['cluster_alg']]
    for result_file in result_files:
        if result_file.endswith(result_format):
            os.rename(result_file,outfile)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','cluster_file',pipeline)

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
