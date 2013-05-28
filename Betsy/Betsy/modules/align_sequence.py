#align_sequence.py
import os
from Betsy import module_utils
from time import strftime,localtime
import subprocess
from genomicode import config


def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    species = parameters['ref']
    if species == 'human':
       ref_file = config.human_ref
    elif species == 'fly':
       ref_file= config.fly_ref
    assert os.path.exists(ref_file),'the ref_file %s does not exist' % ref_file
    outfile = get_outfile(parameters, objects, pipeline)
    command = ['bwa','aln',ref_file,single_object.identifier]
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
    dirname = os.path.dirname(outfile)
    a=os.listdir(dirname)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for align_sequence does not exist' % outfile)
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
    filename = 'align_sequence' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'fastq_file', 'contents,read')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for align_sequence does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    attributes = parameters.values()
    new_object = module_utils.DataObject('sai_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects



