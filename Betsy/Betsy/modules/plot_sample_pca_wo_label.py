#plot_sample_pca_wo_label.py
import os
#import matplotlib.cm as cm
from Betsy import bie3
from Betsy import rulebase
from Betsy import read_label_file
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    module_utils.plot_pca(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for pca_sample_plot fails' % outfile
    )
    out_node = bie3.Data(rulebase.PcaPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = ('Pca_' + original_file + '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='PcaAnalysis')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node
