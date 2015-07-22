#generate_alignment_sam_pair.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    fastq1_node, sai1_node, fastq2_node, sai2_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    species = out_attributes['ref']
    if species == 'hg18':
        ref_file = config.hg18_ref
    elif species == 'hg19':
        ref_file = config.hg19_ref
    elif species == 'dm3':
        ref_file = config.dm3_ref
    elif species == 'mm9':
        ref_file = config.mm9_ref
    else:
        raise ValueError('cannot handle %s' % species)
    assert os.path.exists(ref_file), 'the ref_file %s does not exist' % ref_file
    bwa_BIN = config.bwa
    bwa_module = module_utils.which(bwa_BIN)
    assert bwa_module, 'cannot find the %s' % bwa_BIN
    command = [bwa_BIN, 'sampe', ref_file, sai1_node.identifier,
               sai2_node.identifier, fastq1_node.identifier,
               fastq1_node.identifier]
    f = file(outfile, 'w')
    try:
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=f,
                                   stderr=subprocess.PIPE)
    finally:
        f.close()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for generate_alignment_sam_pair does not exist' %
        outfile)
    out_node = bie3.Data(rulebase.SamFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'generate_alignment_sam' + original_file + '.sam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(
        datatype_name='FastqFile', read="pair1")
    filter2 = module_utils.AntecedentFilter(
        datatype_name='SaiFile', read="pair1")
    filter3 = module_utils.AntecedentFilter(
        datatype_name='FastqFile', read="pair2")
    filter4 = module_utils.AntecedentFilter(
        datatype_name='SaiFile', read="pair2")
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2, filter3,
        filter4)
    return x
