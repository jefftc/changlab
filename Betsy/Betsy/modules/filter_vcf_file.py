#filter_vcf_file.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    bcftools_BIN = config.bcftools
    bcftools_module = module_utils.which(bcftools_BIN)
    assert bcftools_module, 'cannot find the %s' % bcftools_BIN
    vcfutils_BIN = config.vcfutils
    vcfutils_module = module_utils.which(vcfutils_BIN)
    assert bcftools_module, 'cannot find the %s' % bcftools_BIN
    command = [bcftools_BIN, 'view', in_data.identifier, '|', vcfutils_BIN,
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
    out_node = bie3.Data(rulebase.VcfFile, **out_attributes)
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
    filename = 'filter_' + original_file + '.vcf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='VcfFile')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node
