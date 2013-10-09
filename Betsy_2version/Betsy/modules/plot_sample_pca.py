#plot_sample_pca.py
import os
#from Betsy
import module_utils
import shutil
#from Betsy
import read_label_file
import arrayio
import matplotlib.cm as cm
import bie
import rulebase

def run(in_nodes,parameters):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    a,b,c = read_label_file.read(cls_node.attributes['filename'])
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
        module_utils.plot_pca(data_node.attributes['filename'],
                              outfile,opts,legend)
    else:
        module_utils.plot_pca(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for pca_sample_plot fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.PcaPlot,**new_parameters)
    return out_node

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    data_node.attributes['process']
    filename = ('Pca_' + original_file + '_'+
                data_node.attributes['process']+ '.png')
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='PcaAnalysis')
    cls_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='ClassLabelFile')
    return data_node,cls_node




