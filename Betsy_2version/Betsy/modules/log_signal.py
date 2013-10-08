#log_signal.py
import os
import shutil
#from Betsy import module_utils
import module_utils
from genomicode import binreg
import bie
import rulebase

def run(data_node, parameters):
    """log the input file"""
    import arrayio
    import math
    outfile = name_outfile(data_node)
    M = arrayio.read(data_node.attributes['filename'])
    assert not binreg.is_logged_array_data(M), 'the file is logged'
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None:
                if float(M._X[i][j])<1:
                    M._X[i][j] = 1
                M._X[i][j]=math.log(float(M._X[i][j]),2)
    f = file(outfile,'w')
    M_c = arrayio.convert(M,to_format=arrayio.tab_delimited_format)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for log_signal fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_log_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
