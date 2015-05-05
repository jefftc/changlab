#plot_sample_pca.py
import os
from Betsy import module_utils
from Betsy import read_label_file
import matplotlib.cm as cm
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    a, b, c = read_label_file.read(cls_node.identifier)
    if len(a) > 1:
        colors = []
        for i in range(5):
            colors.append(cm.hot(i / 5.0, 1))
            colors.append(cm.autumn(i / 5.0, i))
            colors.append(cm.cool(i / 5.0, i))
            colors.append(cm.jet(i / 5.0, i))
            colors.append(cm.spring(i / 5.0, i))
            colors.append(cm.prism(i / 5.0, i))
            colors.append(cm.summer(i / 5.0, i))
            colors.append(cm.winter(i / 5.0, i))
        opts = [colors[int(i)] for i in b]
        legend = [c[int(i)] for i in b]
        module_utils.plot_pca(data_node.identifier, outfile, opts, legend)
    else:
        module_utils.plot_pca(data_node.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for pca_sample_plot fails' % outfile
    )
    out_node = bie3.Data(rulebase.PcaPlot, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    data_node.attributes['process']
    filename = (
        'Pca_' + original_file + '_' + data_node.attributes['process'] + '.png'
    )
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(
        datatype_name='PcaAnalysis', process=out_attributes["process"])
    filter2 = module_utils.AntecedentFilter(datatype_name='ClassLabelFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x
