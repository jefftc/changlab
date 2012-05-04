#centering.py
import os
import subprocess
import module_utils

def run(parameters,objects,pipeline):
    """mean or median"""
    CLUSTER_BIN = 'cluster'
    center_alg = {'mean':'a','median':'m'}
    try :
        center_parameter = center_alg[parameters['gene_center']]
    except:
        raise ValueError("Centering parameter is not recognized")
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    
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
    assert module_utils.exists_nz(outfile),('the output file %s for centering fails'
                                            %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),('the input file %s for centering does not exist'
                                       %identifier)
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
