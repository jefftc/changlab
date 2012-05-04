#unlog_signal_file.py
import os
import shutil
import module_utils
from genomicode import binreg
def run(parameters,objects,pipeline):
    """unlog the pcl file"""
    import arrayio
    import math
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(identifier)
    assert binreg.is_logged_array_data(M),'the input file\
                               %s should be logged'%identifier
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None :
                M._X[i][j] = 2**float(M._X[i][j])
    f = file(outfile,'w')
    arrayio.pcl_format.write(M,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output\
                        file %s for unlog_signal_file fails'%outfile
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
    assert os.path.exists(identifier),'the input\
            file %s for unlog_signal_file does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

