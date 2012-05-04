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
           'hierarchical':['-m','m','-ag'],'som':["-g",dist,'-s']}
    try:    
        com_parameter = alg[parameters['cluster_alg']]
    except:
        raise ValueError("cluster algorithm is not recognized")
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
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
    result_formats = {'kmeans': 'cdt','pca':'pca_gene.coords.txt',#or 'pc.txt'
           'hierarchical':'nrm','som':'txt'}
    result_format = result_formats[parameters['cluster_alg']]
    for result_file in result_files:
        if result_file.endswith(result_format):
            os.rename(result_file,outfile)
    assert module_utils.exists_nz(outfile),'the output file %s\
               for clustering fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input file %s\
                for clustering does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'cluster_file',parameters,objects,single_object)
    return new_objects
