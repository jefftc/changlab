#call_variants_GATK.py
import os
from Betsy import module_utils
from time import strftime,localtime
import subprocess
from genomicode import config

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    GATK_path = config.Gatk
    GATK_BIN = module_utils.which(GATK_path)
    assert os.path.exists(GATK_path),(
        'cannot find the %s' %GATK_path)
    species = parameters['ref']
    if species == 'hg18':
       ref_file = config.hg18_ref
    elif species == 'hg19':
       ref_file = config.hg19_ref
    elif species == 'dm3':
       ref_file = config.dm3_ref
    elif species == 'mm9':
       ref_file= config.mm9_ref
    assert os.path.exists(ref_file),'the ref file %s does not exsits' %ref_file
    command = ['java','-jar',GATK_path,'-T','UnifiedGenotyper',
              '-R',ref_file,'-I',single_object.identifier,
               '-o',outfile,'-rf', 'BadCigar','-stand_call_conf','50.0',
               '-stand_emit_conf', '10.0']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    process.wait()
    #error_message = process.communicate()[1]
    #print error_message
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for call_variants_GATK does not exist'
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
    filename = 'GATK_' + original_file + '.vcf'
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
        'the input file %s for call_variants_GATK does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'vcf_file', parameters, objects, single_object)
    return new_objects





