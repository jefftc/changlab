#convert_to_tdf_if_not_tdf.py
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
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = None
    tmp_file = None
    try:
        xlrd.open_workbook(identifier)
        tmp_file='tmp.xls'
    except Exception,XLRDError:
        try:
            book = openpyxl.load_workbook(identifier)
            tmp_file = 'tmp.xlsx'
        except Exception,InvalidFileException:
            tmp_file = None
        except (SystemError,MemoryError,KeyError),x:
            raise
        
    if not tmp_file:
        try:
            M = arrayio.choose_format(identifier)
        except Exception,x:
                raise 
        except (SystemError,MemoryError,KeyError),x:
            raise 
    
    if M:
        shutil.copyfile(identifier,outfile)
    elif tmp_file:
        shutil.copyfile(identifier,tmp_file)
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
    assert module_utils.exists_nz(outfile),'the output \
        file %s for convert_to_tdf_if_not_tdf does not exists'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                            single_object,pipeline)
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
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input \
        file %s for convert_to_tdf_if_not_tdf does not exist'%identifier
    return identifier,single_object

