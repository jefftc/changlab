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
    outfile,new_objects = get_outfile(parameters,objects)     
    M = arrayio.read(identifier)
    label_line = ['0'] * M.dim()[1]
    assert parameters['Contents'].startswith('[') and parameters['Contents'].endswith(']')
    class_name=[parameters['Contents'][1:-1]]
    read_label_file.write(outfile,class_name,label_line)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,
                                    objects,'signal_file','Contents,DatasetId')
    
def get_outfile(parameters,objects):
    outfile = os.path.join(os.getcwd(),'class_label_file.cls')
    if 'status' in parameters.keys():
        del parameters['status']
    attributes = parameters.values()
    new_object=rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,
                            objects,'signal_file','Contents,DatasetId')
   
