#center_genes.py
import os
import subprocess
from Betsy import module_utils, bie3, rulebase
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    """mean or median"""
    in_data = antecedents
    CLUSTER_BIN = config.cluster
    cluster = module_utils.which(CLUSTER_BIN)
    assert cluster, 'cannot find the %s' % CLUSTER_BIN
    center_alg = {'mean': 'a', 'median': 'm'}
    try:
        center_parameter = center_alg[out_attributes['gene_center']]
    except:
        raise ValueError("Centering parameter is not recognized")
    
    outfile = name_outfile(in_data, user_options)
    process = subprocess.Popen([CLUSTER_BIN, '-f', in_data.identifier, '-cg',
                                center_parameter, '-u', outfile],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    
    outputfile = outfile + '.nrm'
    os.rename(outputfile, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for centering fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Normalize, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_center_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    new_parameters = out_attributes.copy()
    new_parameters['format'] = 'tdf'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
