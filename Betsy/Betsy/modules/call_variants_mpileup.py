#call_variants_mpileup.py
import os
from Betsy import module_utils
from time import strftime,localtime
import subprocess
from genomicode import config

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    
    ref = config.fly_ref
    ref = '/home/xchen/try_GATK/exampleFiles/exampleFASTA.fasta'
    assert os.path.exists(ref),'the ref file %s does not exsits' %ref
    #command = ['samtools','mpileup','-uf',ref,single_object.identifier,'|',
    #           'bcftools','view', '-bvcg','-','>',outfile]
    command = ['samtools','mpileup','-uf',ref,single_object.identifier]
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
        'the output file %s for call_variants_mpileup does not exist'
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
    filename = 'mpileup_' + original_file + '.bcf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'input_sam_file', 'contents')
    if not single_object:
        single_object = module_utils.find_object(parameters, objects,
                                                 'sam_file',
                                                 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for call_variants_mpileup does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'vcf_file', parameters, objects, single_object)
    return new_objects





