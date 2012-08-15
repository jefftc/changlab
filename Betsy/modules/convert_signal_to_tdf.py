#convert_signal_to_tdf.py
import os
import hash_method
import shutil
import xlrd
import module_utils
import rule_engine
import openpyxl
import arrayio
def run(parameters,objects,pipeline):
    """check an input file is xls or xlsx format"""
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if single_object.identifier.endswith('.gz'):
        newfile = module_utils.gunzip(single_object.identifier)
    else:
        newfile = single_object.identifier
    M = None
    tmp_file = None
    try:
        xlrd.open_workbook(newfile)
        tmp_file='tmp.xls'
    except Exception,XLRDError:
        try:
            book = openpyxl.load_workbook(newfile)
            tmp_file = 'tmp.xlsx'
        except Exception,InvalidFileException:
            tmp_file = None
        except (SystemError,MemoryError,KeyError),x:
            raise
        
    if not tmp_file:
        try:
            M = arrayio.choose_format(newfile)
        except Exception,x:
                raise 
        except (SystemError,MemoryError,KeyError),x:
            raise 
    
    if M:
        shutil.copyfile(newfile,outfile)
    elif tmp_file:
        shutil.copyfile(newfile,tmp_file)
        import Betsy_config
        xls2txt_path = Betsy_config.XLS2TXT
        xls2txt_BIN = module_utils.which(xls2txt_path)
        assert xls2txt_BIN,'cannot find the %s' %xls2txt_path
        f = file(outfile,'w')
        command =['python',xls2txt_BIN, tmp_file]
        process = subprocess.Popen(command,shell=False,
                                    stdout= f ,
                                    stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        os.remove(tmp_file)
        f.close()
    assert module_utils.exists_nz(outfile),(
    'the output file %s for convert_to_tdf_if_not_tdf does not exists'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                            single_object,pipeline,outfile)
    return new_objects


        
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
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    hash_result = make_unique_hash(single_object.identifier,pipeline,parameters)
    filename = 'signal_'+ hash_result + '*' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'input_signal_file','contents')
    assert os.path.exists(single_object.identifier),(
    'the input file %s for convert_to_tdf_if_not_tdf does not exist'
    %single_object.identifier)
    return single_object

