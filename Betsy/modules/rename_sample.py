#rename_sample.py
import os
import shutil
import module_utils
import arrayio
import subprocess
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    rename_file,obj = module_utils.find_object(parameters,objects,
                                          'rename_list_file','contents')
    #module_utils.check_rename_file(rename_file)
    #file_dict = module_utils.read_rename_file(rename_file)
    outfile = get_outfile(parameters,objects,pipeline)
    import Betsy_config
    rename_path = Betsy_config.RENAME
    rename_BIN = module_utils.which(rename_path)
    assert rename_BIN,'cannot find the %s' %rename_path
    command=['python',rename_BIN,identifier,'--relabel_col_ids',
             rename_file+',NewName','-o',outfile]
    #M=arrayio.read(identifier)
    #namelist = M._col_names[M._col_order[0]]
    #for i in range(len(namelist)):
     #   for key in file_dict.keys():
     #       if key in namelist[i]:
     #           namelist[i] = file_dict[key]
     #           break
    #M._col_names[M._col_order[0]]=namelist
    #f = file(outfile,'w')
    #arrayio.pcl_format.write(M,f)
    #f.close()
     
    assert module_utils.exists_nz(outfile),'the output\
        file %s for rename_sample does not exist'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents,preprocess',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(identifier),'the input \
        file %s for median_fill_if_missing does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects 


