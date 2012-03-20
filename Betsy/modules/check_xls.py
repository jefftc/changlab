#check_xls.py
import os
import hash_method
import shutil
import xlrd
import module_utils
import rule_engine
import openpyxl

def run(parameters,objects,pipeline):
    """check an input file is xls or xlsx format"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    try:
        xlrd.open_workbook(identifier)
    except Exception,XLRDError:
        try:
            book = openpyxl.load_workbook(identifier)
        except Exception,InvalidFileException:
            return None
    shutil.copyfile(identifier,outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                            single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    if 'status' in parameters.keys():
        del parameters['status']
    original_file = module_utils.get_inputid(identifier)
    byte_size,md5_checksum,sha1_checksum = hash_method.get_file_checksum(identifier)
    new_parameters = parameters.copy()
    new_parameters['file size']= byte_size
    new_parameters['checksum1']= md5_checksum
    new_parameters['checksum2']= sha1_checksum
    hash_result = hash_method.hash_parameters(original_file,pipeline,**new_parameters)
    return hash_result

    
def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object
