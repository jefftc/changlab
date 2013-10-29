#plot_intensity_boxplot.py
import os
import shutil
import math
from genomicode import mplgraph,jmath
import arrayio
from Betsy import module_utils, bie, rulebase

def run(data_node,parameters,network):
    outfile = name_outfile(data_node)
    M = arrayio.read(data_node.attributes['filename'])
    data = jmath.transpose(M._X)
    tickname=M._col_names['_SAMPLE_NAME']
    fig = mplgraph.boxplot(data,xlabel='Sample Name',ylabel='Signal',
                           title='Signal Intensity',box_label=tickname)
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_intensity_boxplot fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.IntensityPlot,**new_parameters)
    return out_node




def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'intensity_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    
    return data_node
