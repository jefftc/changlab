# select_first_n_genes.py
import os
#from Betsy
import module_utils
from time import strftime,localtime
import bie
import rulebase

def run(data_node,parameters):
    """select a num of genes according to num_features"""
    import arrayio
    outfile = name_outfile(data_node)
    f = file(outfile,'w')
    num_features = int(parameters['num_features'])
    assert num_features>0, 'the num_features should be >0'
    M = arrayio.read(data_node.attributes['filename'])
    M_c = M.matrix(range(0,num_features),None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for select_first_n_genes fails'%outfile)
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
    filename = 'signal_select_n_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):

    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

