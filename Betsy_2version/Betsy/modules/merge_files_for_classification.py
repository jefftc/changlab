#merge_files_for_classification.py

import module_utils
import os
import arrayio
import bie
import rulebase

def run(in_nodes, parameters):
    """merge three signal file to generate a joined signal file"""
    merge_node1, merge_node2 = in_nodes
    assert os.path.exists(merge_node1.attributes['filename']),(
    'the merge_file1 %s in merge_data does not exist'%merge_node1.attributes['filename'])
    assert os.path.exists(merge_node2.attributes['filename']),(
    'the merge_file2 %s in merge_data does not exist'%merge_node2.attributes['filename'])
    outfile = name_outfile(in_nodes)
    file1,file2 = module_utils.convert_to_same_platform(merge_node1.attributes['filename'],
                                                        merge_node2.attributes['filename'])
    f = file(outfile,'w')
    module_utils.merge_two_files(file1,file2,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for merge_files_for_classification fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
    return out_node
    



def name_outfile(in_nodes):
    data_node1,data_node2 = in_nodes
    original_file = module_utils.get_inputid(data_node1.attributes['filename'])
    filename = 'signal_merge' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node1,data_node2 = in_nodes
    identifier = data_node1.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='class0,class1')
    data_node2 = module_utils.get_identifier(network, module_id, data_nodes,
                                           contents='test')
    return data_node1, data_node2
