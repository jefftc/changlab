#plot_illu_biotin_line.py
import os
from Betsy import module_utils, bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    module_utils.plot_line_keywd(in_data.identifier, 'biotin', outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_illu_biotin_line fails' % outfile
    )
    out_node = bie3.Data(rulebase.BiotinPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'biotin_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='ControlFile')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node
