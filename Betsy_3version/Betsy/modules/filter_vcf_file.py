#filter_vcf_file.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(data_node, parameters, user_input, network, num_cores):
    outfile = name_outfile(data_node, user_input)
    bcftools_BIN = config.bcftools
    bcftools_module = module_utils.which(bcftools_BIN)
    assert bcftools_module, 'cannot find the %s' % bcftools_BIN
    vcfutils_BIN = config.vcfutils
    vcfutils_module = module_utils.which(vcfutils_BIN)
    assert bcftools_module, 'cannot find the %s' % bcftools_BIN
    command = [bcftools_BIN, 'view', data_node.identifier, '|', vcfutils_BIN,
               'varFilter', '-D500']
    #command = ['vcfutils.pl','varFilter','-D100',single_object.identifier]
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
        'the output file %s for filter_vcf_file does not exist' % outfile
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
    filename = 'filter_' + original_file + '.vcf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes,
                                            datatype='VcfFile')

    return data_node
