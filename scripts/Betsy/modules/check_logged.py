#check_logged.py
import os
import shutil
import module_utils

def run(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    if is_logged(parameters,objects):
        shutil.copyfile(identifier,outfile)
        module_utils.write_Betsy_parameters_file(parameters,single_object)
        return new_objects
    else:
        return None


def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')

def is_logged(parameters,objects):
    import arrayio
    from genomicode import binreg
    identifier,single_object = get_identifier(parameters,objects)
    M = arrayio.read(identifier)
    return binreg.is_logged_array_data(M)
