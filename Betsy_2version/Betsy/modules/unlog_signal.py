#unlog_signal.py
import os
import shutil
#from Betsy
import module_utils
from genomicode import binreg
import bie
import rulebase

def run(data_node,parameters):
    """unlog the pcl file"""
    import arrayio
    import math
    outfile = name_outfile(data_node)
    M = arrayio.read(data_node.attributes['filename'])
    assert binreg.is_logged_array_data(M),(
        'the input file %s should be logged'%data_node.attributes['filename'])
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None :
                M._X[i][j] = 2**float(M._X[i][j])
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for unlog_signal_file fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_unlog' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):

    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


