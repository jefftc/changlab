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
    outfile,new_objects = get_outfile(parameters,objects)
    try:
        xlrd.open_workbook(identifier)
        shutil.copyfile(identifier,outfile)
        module_utils.write_Betsy_parameters_file(parameters,
                                                 single_object,pipeline)
        return new_objects
    except Exception,XLRDError:
        try:
            book =openpyxl.load_workbook(identifier)
            shutil.copyfile(identifier,outfile)
            module_utils.write_Betsy_parameters_file(parameters,
                                                     single_object,pipeline)
            return new_objects
        except(SystemError,MemoryError,KeyError),x:
            raise
        except Exception,x:
            return None
    except (SystemError,MemoryError,KeyError),x:
            raise 
    except Exception,x:
        return None
    
def make_unique_hash(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    filename = os.path.split(identifier)[-1]
    byte_size,md5_checksum,sha1_checksum = hash_method.get_file_checksum(identifier)
    new_parameters = parameters.copy()
    new_parameters['filename'] = filename
    new_parameters['file size']= byte_size
    new_parameters['checksum1']= md5_checksum
    new_parameters['checksum2']= sha1_checksum
    hash_result = hash_method.hash_parameters(**new_parameters)
    return hash_result

    
def get_outfile(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in filename: 
        original_file = '_'.join(filename.split('_')[:-2])
    else:
        original_file = filename
    hash_string = make_unique_hash(parameters,objects)
    filename = original_file + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    objecttype='signal_file'
    if 'Status' in parameters.keys():
        del parameters['Status']
    attributes = parameters.values()
    new_object = rule_engine.DataObject(objecttype,attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
