#recalibrate_base_quality_score.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    GATK_path = config.Gatk
    assert GATK_path, 'cannot find the %s' % GATK_path
    species = out_attributes['ref']
    if species == 'hg18':
        ref_file = config.hg18_ref
        dbsnp_file = config.hg18_dbsnp
        indels_file = config.hg18_indels
    elif species == 'hg19':
        ref_file = config.hg19_ref
        dbsnp_file = config.hg19_dbsnp
        indels_file = config.hg19_indels
    else:
        raise ValueError('we cannot process %s' % species)
    
    assert os.path.exists(ref_file), 'the ref file %s does not exsits' % ref_file
    assert os.path.exists(
        dbsnp_file), 'the dbsnp file %s does not exsits' % dbsnp_file
    assert os.path.exists(
        indels_file), 'the indels file %s does not exsits' % indels_file
    command = ['java', '-jar', GATK_path, '-T', 'BaseRecalibrator', '-R',
               ref_file, '-I', in_data.identifier, '-knownSites', dbsnp_file,
               '-knownSites', indels_file, '-o', 'recal.bam']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    #error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    process.wait()
    command = ['java', '-jar', GATK_path, '-T', 'PrintReads', '-R', ref_file,
               '-BQSR', 'recal.bam', '-I', in_data.identifier, '-o', outfile]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for base_quality_score_recalibration does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.BamFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'base_quality_score_recalibration_' + original_file + '.bam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='BamFile')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node