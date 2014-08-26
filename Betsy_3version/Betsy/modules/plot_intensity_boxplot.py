#plot_intensity_boxplot.py
import os
import shutil
import math
from genomicode import mplgraph,jmath
import arrayio
from Betsy import module_utils, bie3, rulebase

def run(data_node,parameters,user_input, network):
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    data = jmath.transpose(M._X)
    tickname=M._col_names['_SAMPLE_NAME']
    fig = mplgraph.boxplot(data,xlabel='Sample Name',ylabel='Signal',
                           title='Signal Intensity',box_label=tickname)
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_intensity_boxplot fails'%outfile)
    out_node = bie3.Data(rulebase.IntensityPlot,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'intensity_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    
    return data_node
