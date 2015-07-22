#plot_intensity_boxplot.py
import os
from genomicode import mplgraph, jmath
import arrayio
from Betsy import module_utils, bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    data = jmath.transpose(M._X)
    tickname = M._col_names['_SAMPLE_NAME']
    fig = mplgraph.boxplot(data,
                           xlabel='Sample Name',
                           ylabel='Signal',
                           title='Signal Intensity',
                           box_label=tickname)
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_intensity_boxplot fails' % outfile
    )
    out_node = bie3.Data(rulebase.IntensityPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'intensity_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)

    return data_node
