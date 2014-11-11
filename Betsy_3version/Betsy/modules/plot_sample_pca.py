#plot_sample_pca.py
import os
from Betsy import module_utils
import shutil
from Betsy import read_label_file
import arrayio
import matplotlib.cm as cm
from Betsy import bie3
from Betsy import rulebase

def run(in_nodes,parameters, user_input, network,num_cores):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    a,b,c = read_label_file.read(cls_node.identifier)
    if len(a)>1:
        colors = []
        for i in range(5):
            colors.append(cm.hot(i/5.0,1))
            colors.append(cm.autumn(i/5.0,i))
            colors.append(cm.cool(i/5.0,i))
            colors.append(cm.jet(i/5.0,i))
            colors.append(cm.spring(i/5.0,i))
            colors.append(cm.prism(i/5.0,i))
            colors.append(cm.summer(i/5.0,i))
            colors.append(cm.winter(i/5.0,i))
        opts = [colors[int(i)] for i in b]
        legend = [c[int(i)] for i in b]
        module_utils.plot_pca(data_node.identifier,
                              outfile,opts,legend)
    else:
        module_utils.plot_pca(data_node.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for pca_sample_plot fails'%outfile)
    out_node = bie3.Data(rulebase.PcaPlot,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    data_node.attributes['process']
    filename = ('Pca_' + original_file + '_'+
                data_node.attributes['process']+ '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='PcaAnalysis',
                                            optional_key='process',
                                            optional_value=parameters['process'])
    cls_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='ClassLabelFile')
    return data_node,cls_node




