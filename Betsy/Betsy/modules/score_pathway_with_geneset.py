#score_pathway_with_geneset.py
import os
import subprocess
from Betsy import module_utils,config


def run(parameters, objects, pipeline):
    """analyze geneset"""
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    score_geneset_path = config.SCORE_GENE
    score_geneset_BIN = module_utils.which(score_geneset_path)
    assert score_geneset_BIN,'cannot find the %s' %score_geneset_path
    geneset_object = module_utils.find_object(
        parameters, objects, 'geneset_file', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the geneset_file %s for score_pathway_with_geneset does not exist'
        % single_object.identifier)
    geneset = parameters['geneset']
    allgenes = parameters['allgenes']
    automatch = parameters['automatch']
    command = ['python', score_geneset_BIN, '-o', outfile,
               '--geneset_file', geneset_object.identifier,
               single_object.identifier]
    if allgenes == 'yes_allgenes':
        command.append('--all')
    if automatch == 'yes_automatch':
        command.append('--automatch')
    genesets = geneset.split('/')
    for gene in genesets:
        command.extend(['-g', gene])
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for score_pathway_with_geneset fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'score_geneset_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for score_pathway_with_geneset does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'geneset_analysis', parameters, objects, single_object)
    return new_objects

