#make_class_label_file.py
import module_utils
import os
import read_label_file
import shutil
import rule_engine

def run(parameters,objects,pipeline):
    """generate the class_label_file for signal data"""
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(identifier)
    label_line = ['0'] * M.dim()[1]
    assert parameters['contents'].startswith('[') and parameters['contents'].endswith(']')
    class_name=[parameters['contents'][1:-1]]
    read_label_file.write(outfile,class_name,label_line)
    assert module_utils.exists_nz(outfile),'the \
            output file %s for make_class_label_file fails'%outfile

    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
                    identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    outfile = os.path.join(os.getcwd(),'class_label_file.cls')
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    attributes = parameters.values()
    new_object = rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,
                            objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input \
            file %s for make_class_label_file does not exist'%identifier
    return identifier,single_object
