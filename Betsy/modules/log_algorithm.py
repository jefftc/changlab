#log_algorithm.py

import os
import module_utils

def run(parameters,objects,pipeline):
    """log the input gct or pcl file"""
    import arrayio
    import math
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    f = file(outfile,'w')
    M_format = arrayio.choose_format(identifier)
    M = arrayio.read(identifier)
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None :
                if float(M._X[i][j])<1:
                    M._X[i][j] = 1
                M._X[i][j]=math.log(float(M._X[i][j]),2)
    if M_format.__name__ == 'arrayio.gct_format':
        M_c = arrayio.convert(M,to_format=arrayio.gct_format)
        arrayio.gct_format.write(M_c,f)
    elif M_format.__name__ == 'arrayio.pcl_format':
        M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects


def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','signal_file',pipeline)
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
