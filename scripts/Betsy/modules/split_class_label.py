#split_label_file.py
import read_label_file
import module_utils
import shutil
import rule_engine
import os
import json

def run(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    result,label_line,second_line=read_label_file.read(identifier)
    assert parameters['Contents'].startswith('[') and parameters['Contents'].endswith(']')
    contents=parameters['Contents'][1:-1].split(',')
    content_index=[]
    second_line=[class_name.lower() for class_name in second_line]
    for content in contents:
        try:
            a=second_line.index(content)
        except ValueError:
            return None
        content_index.append(a)
    new_label_line=[]
    for label in label_line:
        newline=[str(i) for i in range(len(content_index))
                 if int(label)==content_index[i]]
        new_label_line.extend(newline)
    read_label_file.write(outfile,contents,new_label_line)
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,
                    objects,'class_label_file','PreContents,PreDatasetid')
    
def get_outfile(parameters,objects):
    outfile = os.path.join(os.getcwd(),'class_label_file.cls')
    attributes = parameters.values()
    new_object=rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,
                                    objects,'class_label_file','PreContents,PreDatasetid')

