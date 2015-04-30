#call_variants_GATK.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(data_node, parameters, user_input, network, num_cores):
    outfile = name_outfile(data_node, user_input)
    GATK_path = config.Gatk
    GATK_BIN = module_utils.which(GATK_path)
    assert os.path.exists(GATK_path), ('cannot find the %s' % GATK_path)
    species = parameters['ref']
    if species == 'hg18':
        ref_file = config.hg18_ref
    elif species == 'hg19':
        ref_file = config.hg19_ref
    elif species == 'dm3':
        ref_file = config.dm3_ref
    elif species == 'mm9':
        ref_file = config.mm9_ref
    assert os.path.exists(ref_file), 'the ref file %s does not exsits' % ref_file
    command = ['java', '-jar', GATK_path, '-T', 'UnifiedGenotyper', '-R',
               ref_file, '-I', data_node.identifier, '-o', outfile, '-rf',
               'BadCigar', '-stand_call_conf', '50.0', '-stand_emit_conf',
               '10.0']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    #process.wait()
    error_message = process.communicate()[1]
    #print error_message
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for call_variants_GATK does not exist' % outfile
    )
    out_node = bie3.Data(rulebase.VcfFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def get_out_attributes(parameters, data_object):
    return parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'GATK_' + original_file + '.vcf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes,
                                            datatype='BamFile')

    return data_node
