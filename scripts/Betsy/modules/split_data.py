#split_data.py

import module_utils
import os
import arrayio
import hash_method
import rule_engine
import shutil
import read_label_file

def run(parameters,objects):
    """extract one signal file to another signal file according to contents"""
    identifier,single_object = get_identifier(parameters,objects)
    class_label_file,obj=module_utils.find_object(parameters,
                            objects,'class_label_file','PreContents,PreDatasetid')
    assert os.path.exists(class_label_file)
    result,label_line,second_line=read_label_file.read(class_label_file)
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
    sample_index=[]
    for i in range(len(label_line)):
        newline=[i for j in range(len(content_index))
                 if int(label_line[i])==content_index[j]]
        sample_index.extend(newline)
    outfile,new_objects = get_outfile(parameters,objects)
    M = arrayio.read(identifier)
    M_c=M.matrix(None,sample_index)
    f_out = file(outfile,'w')
    arrayio.pcl_format.write(M_c,f_out)
    f_out.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','PreContents,PreDatasetid')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','PreContents,PreDatasetid','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','PreContents,PreDatasetid')


