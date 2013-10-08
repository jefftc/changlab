#fill_missing_with_zeros.py
import os
import shutil
#from Betsy
import module_utils
import arrayio
from time import strftime,localtime
import bie
import rulebase

def run(data_node,parameters):
    outfile = name_outfile(data_node)
    assert module_utils.is_missing(data_node.attributes['filename']),'no missing values'
    M = arrayio.read(data_node.attributes['filename'])
    f_out = file(outfile, 'w')
    for i in range(M.dim()[0]):
        for j in range(M.dim()[1]):
            if M._X[i][j] is None:
                M._X[i][j] = '0'
    arrayio.tab_delimited_format.write(M, f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for zero_fill_if_missing fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node



def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def get_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_zero_fill_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


