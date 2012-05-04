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
    class_label_file,obj = module_utils.find_object(parameters,
                        objects,'class_label_file','contents')
    assert os.path.exists(class_label_file)
    identifier,single_object = get_identifier(parameters,objects)
    result,label_line,second_line=read_label_file.read(class_label_file)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(identifier)
    if M.dim()[1]==len(label_line):
        shutil.copyfile(identifier,outfile)
        assert module_utils.exists_nz(outfile),'the output\
                           file %s for pull_out_dataset fails'%outfile
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline)
        return new_objects
    else:
        raise ValueError('the pull_out_dataset fails')
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(parameters,
        objects,'signal_file','[unknown]',pipeline)
    

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,
                                objects,'signal_file','[unknown]')
    assert os.path.exists(identifier),'the input \
            file %s for pull_out_dataset does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
