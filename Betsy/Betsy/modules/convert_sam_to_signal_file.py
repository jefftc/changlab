#convert_sam_to_bam.py
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
    format_type = parameters['format']
    options = None
    if format_type =='sam':
        options = '--sam'
    elif format_type == 'bam':
        options = '--bam'
    if species in ['hg18','hg19']:
        ref_file = config.rna_hum
    elif species == 'dm3':
        ref_file =conf.rna_mouse
    command = ['rsem-calculate-expression',options,single_object.identifier,
               ref_file, '--no-bam-output',species,'-p','8']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    process.wait()
    slice_BIN = config.slice_matrix
    command = ['python',slice_BIN,ref_file+'.genes.results',
               '--select_col_ids','transcript_id,gene_id,TPM']
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
        'the output file %s for convert_sam_to_bam does not exist'
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
    filename = original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'rna_fastq_file', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for convert_fastq_to_signal_file does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'signal_raw', parameters, objects, single_object)
    return new_objects





