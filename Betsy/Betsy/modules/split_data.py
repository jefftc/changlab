#split_data.py

import module_utils
import os
import arrayio
import hash_method
import shutil
import read_label_file

def run(parameters,objects,pipeline):
    """extract one signal file to another signal file according to contents"""
    single_object = get_identifier(parameters,objects)
    class_label_file = module_utils.find_object(parameters,
                            objects,'class_label_file','precontents')
    assert os.path.exists(class_label_file.identifier),(
        'class_label_file %s does not exist'%class_label_file.identifier)
    result,label_line,second_line = read_label_file.read(class_label_file.identifier)
    assert parameters['contents'].startswith('[') and parameters['contents'].endswith(']')
    contents = parameters['contents'][1:-1].split(',')
    content_index = []
    second_line = [class_name.lower() for class_name in second_line]
    for content in contents:
        try:
            a = second_line.index(content)    
        except ValueError:
            raise ValueError('the content value in parameters and cls file does not match')
            return None
        content_index.append(a)
    sample_index = []
    for i in range(len(label_line)):
        newline = [i for j in range(len(content_index))
                 if int(label_line[i]) == content_index[j]]
        sample_index.extend(newline)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    M_c = M.matrix(None,sample_index)
    f_out = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for split_data fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
        parameters,['precontents','status'])
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_split_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(
        parameters,['precontents','status'])
    attributes = parameters.values()
    new_object = module_utils.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','precontents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for split_data does not exist'
        %single_object.identifier)
    return single_object

