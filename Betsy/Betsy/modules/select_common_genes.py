#select_common_genes.py
import sys
import arrayio
import module_utils
import os
from genomicode import matrixlib

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    test_file = module_utils.find_object(
        parameters,objects,'signal_file','testcontents')
    training_file = module_utils.find_object(
        parameters,objects,'signal_file','traincontents')
    assert os.path.exists(training_file.identifier),(
        'the training file %s for select_common_genes does not exists'
        %training_file.identifier)
    assert os.path.exists(test_file.identifier),(
        'the test file %s for select_common_genes does not exists'
        %test_file.identifier)
    file1,file2 = module_utils.convert_to_same_platform(
        training_file.identifier,test_file.identifier)
    training = arrayio.read(file1)
    test = arrayio.read(file2)
    [M_A,M_B] = matrixlib.align_rows(training,test)
    assert M_A.nrow()>0,(
        'there is no common genes betwee %s and %s'
        %(training_file.identifier,test_file.identifier))
    f = file(outfile,'w')
    if parameters['contents'] == parameters['traincontents']:
        arrayio.pcl_format.write(M_A,f)
    elif parameters['contents'] == parameters['testcontents']:
        arrayio.tab_delimited_format.write(M_B,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for select_common_genes fails'%outfile)
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
    filename = 'signal_commom_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for select_common_genes does not exist'
        %single_object.identifier)
    return single_object
def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status',
                'traincontents','testcontents'])
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


        
