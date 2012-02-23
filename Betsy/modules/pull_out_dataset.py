#pull_out_data.py
import module_utils
import os
import arrayio
import hash_method
import rule_engine
import shutil
import read_label_file
def run(parameters,objects,pipeline):
    """pull out the signal file if the class label file is given"""
    class_label_file,obj=module_utils.find_object(parameters,
                        objects,'class_label_file','Contents,DatasetId')
    assert os.path.exists(class_label_file)
    identifier,single_object = get_identifier(parameters,objects)
    result,label_line,second_line=read_label_file.read(class_label_file)
    outfile,new_objects = get_outfile(parameters,objects)
    M=arrayio.read(identifier)
    if M.dim()[1]==len(label_line):
        shutil.copyfile(identifier,outfile)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
        return new_objects
    else:
        return None
    
def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,
                    objects,'signal_file','[unknown],DatasetId')
    
def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,
                    objects,'signal_file','[unknown],DatasetId','signal_file')
    

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,
                                objects,'signal_file','[unknown],DatasetId')
