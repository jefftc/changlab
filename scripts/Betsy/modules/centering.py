#centering.py
import os
import subprocess
import module_utils

def run(parameters,objects):
    """mean or median"""
    CLUSTER_BIN = 'cluster'
    center_alg = {'mean':'a','median':'m'}
    try :
        center_parameter = center_alg[parameters['gene_center']]
    except:
        raise ValueError("Centering parameter is not recognized")
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    process = subprocess.Popen([CLUSTER_BIN,'-f',identifier,
                                '-cg',center_parameter,'-u',outfile],
                                shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    outputfile = outfile + '.nrm'
    os.rename(outputfile,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

    
def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')


