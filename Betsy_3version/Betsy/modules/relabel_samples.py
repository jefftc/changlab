#relabel_samples.py
import os
import shutil
import arrayio
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(in_nodes, parameters, user_input, network, num_cores):
    data_node, rename_node = in_nodes
    outfile = name_outfile(in_nodes, user_input)
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
    out_node = bie3.Data(rulebase._SignalFile_Annotate, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes,
                                            datatype='_SignalFile_Annotate')
    rename_node = module_utils.get_identifier(network, module_id, data_nodes,
                                              user_attributes,
                                              datatype='RenameFile')
    return data_node, rename_node


def name_outfile(in_nodes, user_input):
    data_node, rename_node = in_nodes
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_rename_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, in_nodes):
    data_node, rename_node = in_nodes
    new_parameters = data_node.data.attributes.copy()
    new_parameters['rename_sample'] = 'yes'
    return new_parameters


def make_unique_hash(in_nodes, pipeline, parameters, user_input):
    data_node, rename_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)
