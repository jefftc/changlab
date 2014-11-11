# select_first_n_genes.py
import os
from Betsy import module_utils
from Betsy import bie3
from Betsy import rulebase

def run(data_node,parameters, user_input, network,num_cores):
    """select a num of genes according to num_features"""
    import arrayio
    outfile = name_outfile(data_node,user_input)
    f = file(outfile,'w')
    num_features = 500
    if 'num_features_value' in user_input:
    	num_features = int(user_input['num_features_value'])
    assert num_features>0, 'the num_features should be >0'
    M = arrayio.read(data_node.identifier)
    M_c = M.matrix(range(0,num_features),None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for select_first_n_genes fails'%outfile)
    out_node = bie3.Data(rulebase._SignalFile_Filter,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_select_n_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):

    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

