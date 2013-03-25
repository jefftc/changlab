#add_crossplatform_probeid.py

import os
import shutil
import subprocess
import arrayio
from Betsy import module_utils
from genomicode import jmath, Matrix, arrayplatformlib, config
from time import strftime,localtime

def run(parameters, objects, pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    DATA = arrayio.read(single_object.identifier)
    chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
    platform = parameters['platform']
    assert arrayplatformlib.get_bm_attribute(platform), (
        'the desire platform %s is not recognized by Betsy' % platform)
    if chipname == platform:
        shutil.copyfile(single_object.identifier, outfile)
    else:
        Annot_path = config.annotate_matrix
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN, 'cannot find the %s' % Annot_path
        command = ['python', Annot_BIN, '-f', single_object.identifier,
                   '-o', outfile, "--platform", platform]
        process = subprocess.Popen(command, shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for add_crossplatform_probeid fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(parameters, single_object,
                                             pipeline, outfile,starttime,user,jobname)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    platform = parameters['platform']
    filename = 'signal_' + platform + '_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for add_crossplatform_probeid'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'signal_file', parameters, objects, single_object)
    return new_objects
