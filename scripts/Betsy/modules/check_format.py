#check_format.py
import shutil
import module_utils
import arrayio
def run(parameters,objects):
    """check if not_xls format is res"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    M=arrayio.choose_format(identifier)
    if parameters['format']==M.__name__[8:-7]:
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


