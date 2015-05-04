#fix_header_GATK.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    AddOrReplaceReadGroups_path = config.AddOrReplaceReadGroups
    assert os.path.exists(AddOrReplaceReadGroups_path), (
        'cannot find the %s' % AddOrReplaceReadGroups_path
    )
    command = ['java', '-Xmx5g', '-jar', AddOrReplaceReadGroups_path,
               'I=' + in_data.identifier, 'O=' + outfile, 'PL=illumina',
               'ID=Group1', 'LB=Al_chrom3', 'PU=Al_chrom3', 'SM=Al_chrom3',
               'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=LENIENT']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for fix_header_GATK does not exist' % outfile
    )
    out_node = bie3.Data(rulebase.BamFile, **out_attributes)
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
    filename = 'fix_header_' + original_file + '.bam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='BamFile')

    return data_node
