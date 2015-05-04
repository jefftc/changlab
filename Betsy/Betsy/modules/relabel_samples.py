#relabel_samples.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, rename_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    rename_path = config.slice_matrix
    rename_BIN = module_utils.which(rename_path)
    assert rename_BIN, 'cannot find the %s' % rename_path
    command = ['python', rename_BIN, data_node.identifier, '--relabel_col_ids',
               rename_node.identifier + ',NewName']
    f = file(outfile, 'w')
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=f,
                               stderr=subprocess.PIPE)
    f.close()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)

    assert module_utils.exists_nz(outfile), (
        'the output file %s for relabel_samples does not exist' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Annotate, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='_SignalFile_Annotate')
    rename_node = module_utils.get_identifier(network, module_id, pool,
                                              user_attributes,
                                              datatype='RenameFile')
    return data_node, rename_node


def name_outfile(antecedents, user_options):
    data_node, rename_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_rename_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    data_node, rename_node = antecedents
    new_parameters = data_node.data.attributes.copy()
    new_parameters['rename_sample'] = 'yes'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, rename_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
