#dwd.py

import os
import module_utils
import shutil
import read_label_file
from genomicode import dwdnorm
import arrayio
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),'cannot find label_file %s'%label_file.identifier
    M = arrayio.read(single_object.identifier)
    result,label_line,second_line=read_label_file.read(label_file.identifier)
    assert len(result) == 2, 'for dwd,there should be only 2 classes'
    assert [i in ['0','1'] for i in label_line] == [True]*len(label_line),(
        'the label of class shoul be 0 and 1')
    y = [i.replace('0','-1') for i in label_line]
    M_y = dwdnorm.normalize(M,y)
  
    f = file(outfile,'w')
    arrayio.pcl_format.write(M_y,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for dwd fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects


def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_dwd_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for dwd does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
