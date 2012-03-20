#convert_xls_pcl.py
import os
import module_utils
import xlrd
import openpyxl
import shutil
import subprocess

def run(parameters,objects,pipeline):
    """convert xls or xlsx signal file to pcl format"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    #convert it to tdf format
    try:
        book = xlrd.open_workbook(identifier)
        tmp_file = 'tmp.xls'
    except Exception, XLRDError:
        tmp_file = 'tmp.xlsx'
    shutil.copyfile(identifier,tmp_file)
    import Betsy_config
    xls2txt_BIN = Betsy_config.XLS2TXT
    f = file(outfile,'w')
    command =['python',xls2txt_BIN, tmp_file]
    process = subprocess.Popen(command,shell=False,
                                stdout= f ,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    os.remove(tmp_file)
    #convert to pcl format and rewrite 
    import arrayio
    M = arrayio.read(outfile)
    M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
    f = file(outfile,'w')
    arrayio.pcl_format.write(M_c,f)
    f.close()
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)
    
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
