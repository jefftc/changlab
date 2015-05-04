#plot_illu_hyb_bar.py
import os
from genomicode import mplgraph
import numpy
import math
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    plot_hyb_bar(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for plot_hyb_bar fails' % outfile
    )
    out_node = bie3.Data(rulebase.Hyb_barPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'hyb_bar_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='ControlFile')

    return data_node


def plot_hyb_bar(filename, outfile):
    high = ['ILMN_2038770', 'ILMN_2038769']
    med = ['ILMN_2038768', 'ILMN_2038771']
    low = ['ILMN_1343050', 'ILMN_1343052']
    high_data = []
    med_data = []
    low_data = []
    import arrayio
    M = arrayio.read(filename)
    header = M.row_names()
    for i in range(M.dim()[0]):
        if not M.row_names(header[1])[i] == 'cy3_hyb':
            continue
        if M.row_names(header[0])[i] in high:
            high_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in med:
            med_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in low:
            low_data.extend(M.slice()[i])
    mean = [numpy.mean(high_data), numpy.mean(med_data), numpy.mean(low_data)]
    flag = [math.isnan(i) for i in mean]
    assert True not in flag, 'input is not a control file'
    std = [numpy.std(high_data), numpy.std(med_data), numpy.std(low_data)]
    fig = mplgraph.barplot(mean, std,
                           ylabel='Signal',
                           box_label=['high', 'med', 'low'])
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile), 'the plot_illu_hyb_bar.py fails'
