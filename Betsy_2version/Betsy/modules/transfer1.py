#transfer1.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie, rulebase

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    parameters = get_out_attributes(parameters,data_node)
    shutil.copyfile(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for log_signal fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile1,**new_parameters)
    return out_node



def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_file1_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
