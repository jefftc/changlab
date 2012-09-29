#make_class_label_file.py
from Betsy import module_utils
import os
from Betsy import read_label_file
import shutil

def run(parameters,objects,pipeline):
    """generate the class_label_file for signal data"""
    import arrayio
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    label_line = ['0'] * M.dim()[1]
    assert parameters['contents'].startswith('[') and parameters['contents'].endswith(']')
    class_name=[parameters['contents'][1:-1]]
    read_label_file.write(outfile,class_name,label_line)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for make_class_label_file fails'%outfile)

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
    filename = 'class_label_file_'+original_file+'.cls'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'class_label_file',parameters,objects,single_object)
    return new_objects

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(parameters,
                            objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for make_class_label_file does not exist'
        %single_object.identifier)
    return single_object
