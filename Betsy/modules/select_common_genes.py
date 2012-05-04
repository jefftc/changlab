#select_common_genes.py
import sys
import arrayio
import module_utils
import os
from genomicode import binreg
import rule_engine

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    test_file,test_obj = module_utils.find_object(
        parameters,objects,'signal_file','testcontents')
    training_file,train_obj = module_utils.find_object(
        parameters,objects,'signal_file','traincontents')
    assert os.path.exists(training_file),'the training file %s\
                     for select_common_genes does not exists' %training_file
    assert os.path.exists(test_file),'the test file %s\
                    for select_common_genes does not exists' %test_file
    training = arrayio.read(training_file)
    test = arrayio.read(test_file)
    [M_A,M_B] = binreg.align_rows(training,test)
    assert M_A.nrow()>0,'there is no common genes betwee %s \
                             and %s'%(training_file,test_file)
    f = file(outfile,'w')
    if parameters['contents'] == parameters['traincontents']:
        arrayio.pcl_format.write(M_A,f)
    elif parameters['contents'] == parameters['testcontents']:
        arrayio.pcl_format.write(M_B,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                           for select_common_genes fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                         parameters,single_object,pipeline)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
                            identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier), 'the input file %s \
        for select_common_genes does not exist' %identifier
    return identifier,single_object
def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status',
                'traincontents','testcontents'])
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


        
