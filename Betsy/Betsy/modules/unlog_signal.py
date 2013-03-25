#unlog_signal.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    """unlog the pcl file"""
    starttime = strftime(module_utils.FMT, localtime())
    import arrayio
    import math
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    assert binreg.is_logged_array_data(M),(
        'the input file %s should be logged'%identifier)
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None :
                M._X[i][j] = 2**float(M._X[i][j])
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for unlog_signal_file fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline,outfile,starttime,user,jobname)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_unlog' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile
    
def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for unlog_signal_file does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'signal_file', parameters, objects, single_object)
    return new_objects

