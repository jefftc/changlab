#check_not_xls.py
import os
import hash_method
import shutil
import xlrd
import rule_engine
import module_utils
def run(parameters,objects,pipeline):
    """check an input file is not xls format"""
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    try:
        M = arrayio.choose_format(identifier)
        if M:
            shutil.copyfile(identifier,outfile)
            new_objects = get_newobjects(parameters,objects,pipeline)
            module_utils.write_Betsy_parameters_file(
                parameters,single_object,pipeline)
            return new_objects
    except (SystemError,MemoryError,KeyError),x:
        raise 
    except Exception,x:
          return None
    
        
    
def make_unique_hash(identifier,pipeline,parameters):
    if 'status' in parameters.keys():
        del parameters['status']
    original_file = module_utils.get_inputid(identifier)
    byte_size,md5_checksum,sha1_checksum = hash_method.get_file_checksum(
                                                              identifier)
    new_parameters = parameters.copy()
    new_parameters['file size']= byte_size
    new_parameters['checksum1']= md5_checksum
    new_parameters['checksum2']= sha1_checksum
    hash_result = hash_method.hash_parameters(
        original_file,pipeline,**new_parameters)
    return hash_result

def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
