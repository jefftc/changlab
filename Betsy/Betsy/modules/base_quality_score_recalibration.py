#base_quality_score_recalibration.py
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
    assert GATK_path,'cannot find the %s' %GATK_path
    species = parameters['ref']
    if species == 'hg18':
       ref_file = config.hg18_ref
       dbsnp_file = config.hg18_dbsnp
       indels_file = config.hg18_indels
    elif species == 'hg19':
       ref_file = config.hg19_ref
       dbsnp_file = config.hg19_dbsnp
       indels_file = config.hg19_indels
    assert os.path.exists(ref_file),'the ref file %s does not exsits' %ref_file
    assert os.path.exists(dbsnp_file),'the dbsnp file %s does not exsits' %dbsnp_file
    assert os.path.exists(indels_file),'the indels file %s does not exsits' %indels_file
    command = ['java','-jar',GATK_path,'-T','BaseRecalibrator','-R',ref_file,
               '-I',single_object.identifier,'-knownSites',dbsnp_file,
               '-knownSites',indels_file,
               '-o','recal.bam']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    #error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    process.wait()
    command = ['java','-jar',GATK_path,'-T','PrintReads','-R',ref_file,
               '-BQSR','recal.bam','-I',single_object.identifier,
               '-o',outfile]
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    process.wait()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for base_quality_score_recalibration does not exist'
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
    filename = 'base_quality_score_recalibration' + original_file + '.bam'
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
        'the input file %s for base_quality_score_recalibration does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'sam_file', parameters, objects, single_object)
    return new_objects





