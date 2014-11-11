#convert_impute_merge.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
import arrayio
from Betsy import bie3, rulebase

def run(data_node,parameters, user_input,network,num_cores):
    outfile = name_outfile(data_node,user_input)
    parameters = get_out_attributes(parameters,data_node)
    shutil.copyfile(data_node.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for convert_postprocess_impute fails'%outfile)
    out_node = bie3.Data(rulebase._SignalFile_Merge,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object



def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_file1_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,pool,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            pool,user_attributes)
    return data_node
