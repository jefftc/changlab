#analyze_phenotype_for_EMT.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cel_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    out_attributes = get_out_attributes(out_attributes, data_node)
    phenotype_BIN = config.analyze_phenotype
    assert os.path.exists(phenotype_BIN)
    assert "geneset_value" in user_options, 'no geneset are provided'
    if not os.path.exists(outfile):
        os.mkdir(outfile)

    command = ['python', phenotype_BIN, '--phenotype', 'EMT',
               '--ignore_samples', 'shCDH1,1', '--gene',
               user_options['geneset_value'], '-o', outfile + '/EMT',
               data_node.identifier, cel_node.identifier]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    x = process.communicate()
    error_message = x[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for analyze_phenotype fails' % outfile
    )
    out_node = bie3.Data(rulebase.EMTAnalysis, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    data_node, cel_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'EMT_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cel_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='SignalFile')
    cel_node = module_utils.get_identifier(network, module_id, pool,
                                           user_attributes,
                                           datatype='CellTypeFile')
    return data_node, cel_node
