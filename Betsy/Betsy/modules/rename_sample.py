#rename_sample.py
import os
import shutil
from Betsy import module_utils
import arrayio
import subprocess
from Betsy import config
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    rename_file = module_utils.find_object(parameters,objects,
                                          'rename_list_file','contents')
    outfile = get_outfile(parameters,objects,pipeline)
    rename_path = config.RENAME
    rename_BIN = module_utils.which(rename_path)
    assert rename_BIN,'cannot find the %s' %rename_path
    command=['python',rename_BIN,single_object.identifier,'--relabel_col_ids',
             rename_file.identifier+',NewName']
    f=file(outfile,'w')
    process = subprocess.Popen(command,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
    f.close()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
     
    assert module_utils.exists_nz(outfile),(
        'the output file %s for rename_sample does not exist'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_rename_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for median_fill_if_missing does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects 


