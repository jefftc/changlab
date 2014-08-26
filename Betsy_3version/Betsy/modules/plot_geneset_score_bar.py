#plot_geneset_score_bar.py

import os
from Betsy import module_utils, bie3, rulebase
import shutil
from genomicode import mplgraph, filelib, jmath


def run(data_node, parameters, user_input, network):
    outfile = name_outfile(data_node,user_input)
    matrix = [x for x in filelib.read_cols(data_node.identifier)]
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
    out_node = bie3.Data(rulebase.GenesetPlot,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'geneset_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node

