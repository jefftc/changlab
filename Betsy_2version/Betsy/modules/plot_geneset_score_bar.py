#plot_geneset_score_bar.py

import os
from Betsy import module_utils, bie, rulebase
import shutil
from genomicode import mplgraph, filelib, jmath


def run(data_node, parameters, network):
    outfile = name_outfile(data_node)
    matrix = [x for x in filelib.read_cols(data_node.attributes['filename'])]
    matrix = [x[1:] for x in matrix]
    matrix = jmath.transpose(matrix)
    sample = matrix[0][1:]
    data = matrix[1:]
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for one_data in data:
        value = one_data[1:]
        value = [float(i) for i in value]
        pair = [(value[i],sample[i]) for i in range(len(value))]
        pair.sort()
        gene_value = [i[0] for i in pair]
        label = [i[1] for i in pair]
        ylabel=one_data[0]
        from genomicode import mplgraph
        fig=mplgraph.barplot(gene_value,box_label=label,xtick_rotation=90,
                             xlabel='sample',ylabel=ylabel)
        output = os.path.join(outfile,ylabel)
        fig.savefig(output+'.png')
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_geneset_score_bar fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.GenesetPlot,**new_parameters)
    return out_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'geneset_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

