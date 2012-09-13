#merge_data.py

import module_utils
import os
import arrayio
import hash_method

def run(parameters,objects,pipeline):
    """merge two signal file to generate a joined signal file"""
    merge_file1 = module_utils.find_object(parameters,objects,'signal_file','merge1,format')
    merge_file2 = module_utils.find_object(parameters,objects,'signal_file','merge2,format')
    assert os.path.exists(merge_file1.identifier),'the merge_file1 %s in merge_data does not exist'%merge_file1
    assert os.path.exists(merge_file2.identifier),'the merge_file2 %s in merge_data does not exist'%merge_file2
    outfile = get_outfile(parameters,objects,pipeline)
    file1,file2 = module_utils.convert_to_same_platform(merge_file1.identifier,merge_file2.identifier)
    f = file(outfile,'w')
    module_utils.merge_two_files(file1,file2,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for merge_data fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,[merge_file1,merge_file2],pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
        parameters,['merge1','merge2','status'])
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)


def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_merge_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','merge1,format')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for merge_data does not exist'%identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(
        parameters,['merge1','merge2','status'])
    attributes = parameters.values()
    new_object = module_utils.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
