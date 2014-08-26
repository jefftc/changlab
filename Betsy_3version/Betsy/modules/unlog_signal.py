#unlog_signal.py
import os
import shutil
from genomicode import binreg
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, user_input,network):
    """unlog the pcl file"""
    import arrayio
    import math
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    assert binreg.is_logged_array_data(M),(
        'the input file %s should be logged'%data_node.identitifer)
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None :
                M._X[i][j] = 2**float(M._X[i][j])
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for unlog_signal_file fails' % outfile)
    out_node = bie3.Data(rulebase.SignalFile_Filter,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_unlog' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):

    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


