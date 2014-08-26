#normalize_samples_with_quantile.py

import os
from genomicode import quantnorm
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, user_input,network):
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    Y = quantnorm.normalize(M)
    f = file(outfile,'w')
    Y_c = arrayio.convert(Y,to_format = arrayio.pcl_format)
    arrayio.tab_delimited_format.write(Y_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for quantile fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile_Merge,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object




def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_quantile_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
    
def get_out_attributes(parameters,data_node):
    new_parameters = data_node.data.attributes.copy()
    new_parameters['quantile_norm']='yes'
    return new_parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)
