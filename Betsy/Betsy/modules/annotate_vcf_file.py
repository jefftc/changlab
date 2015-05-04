#annotate_vcf_file.py
from genomicode import arrayannot, arrayplatformlib, config
from Betsy import module_utils, bie3, rulebase
import os
import subprocess


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    species = out_attributes['ref']
    annotate_BIN = config.annotate_vcf
    command = ['python', annotate_BIN, in_data.identifier, '-o', outfile,
               '-species', species]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annot_vcf_file fails' % outfile
    )
    out_node = bie3.Data(rulebase.VcfFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'vcf_annot_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='VcfFile')

    return data_node
