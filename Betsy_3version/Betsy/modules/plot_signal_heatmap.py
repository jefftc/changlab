#plot_signal_heatmap.py
import os
import subprocess
from genomicode import config, graphlib
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, user_input, network):
    """generate a heatmap of input file"""
    outfile =  name_outfile(data_node,user_input)
    Heatmap_path = config.arrayplot
    Heatmap_BIN = module_utils.which(Heatmap_path)
    assert Heatmap_BIN,'cannot find the %s' %Heatmap_path

    command = ['python', Heatmap_BIN,data_node.identifier,
               '-o',outfile,"--label_arrays",
               "--label_genes",'--no_autoscale']#"--grid"
    if 'color' in parameters.keys():
        color=['--color' , parameters['color'].replace('_','-')]
        command.extend(color)
    M = arrayio.read(data_node.identifier)
    nrow = M.nrow()
    ncol = M.ncol()
    ratio = float(nrow)/ncol
    max_box_height = None
    max_box_width = None
    if 'hm_width_value' in user_input:
        max_box_width = user_input['hm_width_value']
    if 'hm_height_value' in user_input:
         max_box_height = user_input['hm_height_value']
    if ratio >= 4:
        x,y=graphlib.find_tall_heatmap_size(nrow,ncol,
                                            max_box_height=max_box_height,
                                            max_box_width=max_box_width,
                                            max_megapixels=128)
    else:
        x,y=graphlib.find_wide_heatmap_size(nrow,ncol,
                                            max_box_height=max_box_height,
                                            max_box_width=max_box_width,
                                            max_megapixels=128)
    command.extend(['-x',str(x),'-y',str(y)])
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_signal_heatmap fails' %outfile)
    out_node = bie3.Data(rulebase.Heatmap,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'heatmap_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
    
