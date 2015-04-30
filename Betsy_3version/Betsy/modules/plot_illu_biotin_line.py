#plot_illu_biotin_line.py
import os
from Betsy import module_utils, bie3, rulebase
import shutil


def run(data_node, parameters, user_input, network, num_cores):
    outfile = name_outfile(data_node, user_input)
    module_utils.plot_line_keywd(data_node.identifier, 'biotin', outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_illu_biotin_line fails' % outfile
    )
    out_node = bie3.Data(rulebase.BiotinPlot, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'biotin_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_node):
    return parameters


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes,
                                            datatype='ControlFile')

    return data_node
