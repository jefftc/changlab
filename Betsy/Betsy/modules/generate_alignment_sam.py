#generate_alignment_sam.py
import os
from Betsy import module_utils
from time import strftime,localtime
import subprocess
from genomicode import config


def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    species = parameters['ref']
    if species == 'hg18':
       ref_file = config.hg18_ref
    elif species == 'hg19':
       ref_file = config.hg19_ref
    elif species == 'dm3':
       ref_file = config.dm3_ref
    elif species == 'mm9':
       ref_file= config.mm9_ref
    assert os.path.exists(ref_file),'the ref_file %s does not exist' % ref_file
    read = parameters['read']
    if read == 'single':
        fast_file =  module_utils.find_object(
            parameters, objects, 'fastq_file', 'contents,read')
        command = ['bwa','samse',ref_file, single_object.identifier,
                   fast_file.identifier]
    elif read == 'pair':
        newparameters = parameters.copy()
        newparameters['read']='pair1' 
        first_fa = module_utils.find_object(
            newparameters, objects, 'fastq_file', 'contents,read')
        newparameters['read']='pair2' 
        second_fa = module_utils.find_object(
            newparameters, objects, 'fastq_file', 'contents,read')
        second_sai = module_utils.find_object(
            newparameters, objects, 'sai_file', 'contents,read')
        command = ['bwa','sampe',ref_file,single_object.identifier, second_sai.identifier,
                   first_fa.identifier,second_fa.identifier]
    f=file(outfile,'w')
    try:
        process=subprocess.Popen(command,shell=False,
                             stdout=f,
                             stderr=subprocess.PIPE)
    finally:
        f.close()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for generate_alignment_sam does not exist'
        % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile,starttime,user,jobname)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'generate_alignment_sam' + original_file+ '.sam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    read = parameters['read']
    newparameters = parameters.copy()
    if read == 'pair':  
        newparameters['read']='pair1' 
    single_object = module_utils.find_object(
        newparameters, objects, 'sai_file', 'contents,read')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for generate_alignment_sam does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'sam_file', parameters, objects, single_object)
    return new_objects


