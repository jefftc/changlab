#zero_fill.py
import os
import module_utils

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    import arrayio
    f_out = file(outfile,'w')
    M = arrayio.read(identifier)
    for i in range(M.dim()[0]):
       for j in range(M.dim()[1]):
           if M._X[i][j] == None:
                    M._X[i][j] = '0'
    arrayio.pcl_format.write(M,f_out)    
    f_out.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects


def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
