#split_label_file.py
import read_label_file
import module_utils
import shutil
import rule_engine
import os
import json

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    
    result,label_line,second_line = read_label_file.read(identifier)
    assert parameters['contents'].startswith('[') and parameters['contents'].endswith(']')
    contents = parameters['contents'][1:-1].split(',')
    content_index = []
    second_line = [class_name.lower() for class_name in second_line]
    for content in contents:
        try:
            a = second_line.index(content)
        except ValueError:
            return None
        content_index.append(a)
    new_label_line=[]
    for label in label_line:
        newline=[str(i) for i in range(len(content_index))
                 if int(label)==content_index[i]]
        new_label_line.extend(newline)
    read_label_file.write(outfile,second_line,new_label_line)
    assert module_utils.exists_nz(outfile),'the output\
                                file %s for split_class_label fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
             parameters,['precontents'])
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'class_label_file','precontents',pipeline)

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(
             parameters,['precontents'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,
                    objects,'class_label_file','precontents')
    assert os.path.exists(identifier),'the input\
                file %s for split_class_label does not exist'%identifier
    return identifier,single_object
