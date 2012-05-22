#merge_data.py

import module_utils
import os
import arrayio
import hash_method
import rule_engine


def run(parameters,objects,pipeline):
    """merge two signal file to generate a joined signal file"""
    merge_file1,obj1 = module_utils.find_object(parameters,objects,'signal_file','merge1')
    merge_file2,obj2 = module_utils.find_object(parameters,objects,'signal_file','merge2')
    assert os.path.exists(merge_file1),'the merge_file1 %s in merge_data does not exist'%merge_file1
    assert os.path.exists(merge_file2),'the merge_file2 %s in merge_data does not exist'%merge_file2
    outfile = get_outfile(parameters,objects,pipeline)
    f = file(outfile,'w')
    module_utils.merge_two_files(merge_file1,merge_file2,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s for merge_data fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,[obj1,obj2],pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
        parameters,['merge1','merge2','status'])
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)


def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    parameters = module_utils.renew_parameters(
        parameters,['merge1','merge2','status'])
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','merge1')
    assert os.path.exists(identifier),'the input file %s for merge_data does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(
        parameters,['merge1','merge2','status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
