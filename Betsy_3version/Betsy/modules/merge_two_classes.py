#merge_two_class.py

import os
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(in_nodes, parameters, user_input, network):
    """merge three signal file to generate a joined signal file"""
    merge_node1, merge_node2 = in_nodes
    assert os.path.exists(merge_node1.identifier),(
    'the merge_file1 %s in merge_data does not exist'%merge_node1.identifier)
    assert os.path.exists(merge_node2.identifier),(
    'the merge_file2 %s in merge_data does not exist'%merge_node2.identifier)
    outfile = name_outfile(in_nodes,user_input)
    file1,file2 = module_utils.convert_to_same_platform(merge_node1.identifier,
                                                        merge_node2.identifier)
    f = file(outfile,'w')
    module_utils.merge_two_files(file1,file2,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for merge_data fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile_Merge,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def name_outfile(in_nodes,user_input):
    data_node1,data_node2 = in_nodes
    original_file = module_utils.get_inputid(data_node1.identifier)
    filename = 'signal_merge' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node1,data_node2 = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                             contents='class0')
    data_node2 = module_utils.get_identifier(network, module_id, data_nodes,
                                             user_attributes,
                                           contents='class1')
    return data_node1, data_node2
