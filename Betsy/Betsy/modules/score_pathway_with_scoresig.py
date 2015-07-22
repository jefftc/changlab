#score_pathway_with_scoresig.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    scoresig_path = config.scoresig
    scoresig_BIN = module_utils.which(scoresig_path)
    assert scoresig_BIN, 'cannot find the %s' % scoresig_path
    command = ['python', scoresig_BIN, '-r', in_data.identifier, '-m',
               in_data.identifier, '-j', '20', '-o', outfile]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for run_scoresig does not exists' % outfile
    )
    out_node = bie3.Data(rulebase.SignatureScore, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='SignalFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return x


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signature_score' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
