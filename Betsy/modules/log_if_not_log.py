#log_if_not_log.py
import os
import shutil
import module_utils
from genomicode import binreg
def run(parameters,objects,pipeline):
    """log the input gct or pcl file"""
    import arrayio
    import math
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    if binreg.is_logged_array_data(M):
        shutil.copyfile(single_object.identifier,outfile)
    else:
        for i in range(len(M._X)):
            for j in range(len(M._X[i])):
                if M._X[i][j] is not None :
                    if float(M._X[i][j])<1:
                        M._X[i][j] = 1
                    M._X[i][j]=math.log(float(M._X[i][j]),2)
        f = file(outfile,'w')
        M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
        f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for log_if_not_log fails'%outfile)

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
    filename = 'signal_log_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for log_if_not_log does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

