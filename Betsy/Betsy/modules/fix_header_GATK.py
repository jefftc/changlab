#fix_header_GATK.py
import os
from Betsy import module_utils
from time import strftime,localtime
import subprocess
from genomicode import config

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    AddOrReplaceReadGroups_path = config.AddOrReplaceReadGroups
    assert os.path.exists(AddOrReplaceReadGroups_path),('cannot find the %s'
                                        %AddOrReplaceReadGroups_path)
    command = ['java','-Xmx5g','-jar',AddOrReplaceReadGroups_path,
               'I='+single_object.identifier,
               'O='+outfile,'PL=illumina',
               'ID=Group1','LB=Al_chrom3', 'PU=Al_chrom3',
               'SM=Al_chrom3','CREATE_INDEX=true',
               'VALIDATION_STRINGENCY=LENIENT']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    
    error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for fix_header_GATK does not exist'
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
    filename = 'fix_header_' + original_file + '.bam'
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
        'the input file %s for fix_header_GATK does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'sam_file', parameters, objects, single_object)
    return new_objects





