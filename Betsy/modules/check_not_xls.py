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
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    try:
        arrayio.choose_format(identifier)
        shutil.copyfile(identifier,outfile)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
        return new_objects
    except (SystemError,MemoryError,KeyError),x:
        raise 
    except Exception,x:
          return None
    
        
    
def make_unique_hash(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in filename: 
        original_file = '_'.join(filename.split('_')[:-2])
    else:
        original_file = filename
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
    filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in filename: 
        original_file = '_'.join(filename.split('_')[:-2])
    else:
        original_file = filename
    hash_string = make_unique_hash(parameters,objects,pipeline)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    objecttype = 'signal_file'
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

