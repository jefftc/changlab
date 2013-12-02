#normalize_samples_with_quantile.py

import os
from genomicode import quantnorm
import arrayio
from Betsy import bie 
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, network):
    outfile = name_outfile(data_node)
    M = arrayio.read(data_node.attributes['filename'])
    Y = quantnorm.normalize(M)
    f = file(outfile,'w')
    Y_c = arrayio.convert(Y,to_format = arrayio.pcl_format)
    arrayio.tab_delimited_format.write(Y_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for quantile fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_quantile_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
