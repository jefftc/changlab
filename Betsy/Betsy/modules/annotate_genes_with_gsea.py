#annotate_genes_with_gsea.py
import subprocess
import shutil
import os
import arrayio
from genomicode import arrayplatformlib, config
from Betsy import module_utils, read_label_file
from time import strftime,localtime

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    label_file = module_utils.find_object(
        parameters, objects, 'class_label_file', 'contents')
    assert os.path.exists(label_file.identifier), (
        'cannot find label_file %s for gsea' % label_file.identifier)
    gsea_path = config.gsea
    gsea_module = module_utils.which(gsea_path)
    assert gsea_module, 'cannot find the %s' % gsea_path
    M = arrayio.read(single_object.identifier)
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    chipname = x[0][1]
    platform = chipname + '.chip'
    download_directory = os.path.join(os.getcwd(),'gsea_result')
    command = [gsea_module, single_object.identifier, '--cls_file',
               label_file.identifier, '--platform', platform,
               download_directory]
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(download_directory), (
        'there is no output directory for GSEA')
    shutil.copytree(download_directory, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annotate_genes_with_gsea fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, [single_object,label_file], pipeline,
        outfile,starttime,user,jobname)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the train file %s for annotate_genes_with_gsea does not exist'
        % single_object.identifier)
    return single_object


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'gsea_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'signature_analysis', parameters, objects, single_object)
    return new_objects
