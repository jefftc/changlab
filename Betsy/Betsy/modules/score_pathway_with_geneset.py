#score_pathway_with_geneset.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    """analyze geneset"""
    data_node, geneset_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    score_geneset_path = config.score_geneset
    score_geneset_BIN = module_utils.which(score_geneset_path)
    assert score_geneset_BIN, 'cannot find the %s' % score_geneset_path
    geneset = user_options['geneset_value']
    assert geneset, 'please select geneset to score pathway'
    automatch = out_attributes['automatch']
    command = ['python', score_geneset_BIN, '-o', outfile, '--geneset_file',
               geneset_node.identifier, data_node.identifier]
    if automatch == 'yes':
        command.append('--automatch')
    genesets = geneset.split('/')
    for gene in genesets:
        command.extend(['-g', gene])
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for score_pathway_with_geneset fails' % outfile
    )
    out_node = bie3.Data(rulebase.GenesetAnalysis, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='SignalFile')
    filter2 = module_utils.AntecedentFilter(datatype_name='GenesetFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'score_geneset_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
