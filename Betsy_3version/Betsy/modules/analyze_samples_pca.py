#analyze_samples_pca.py
import os
from genomicode import pcalib
import arrayio
from Betsy import module_utils
from Betsy import bie3
from Betsy import rulebase

def run(data_node,parameters,user_input,network):
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    X = M._X
    N = 500
    if 'pca_gene_num' in user_input:
        N = int(user_input['pca_gene_num'])
    N = min(N,M.nrow())
    index = pcalib.select_genes_var(X,N)
    M_new = M.matrix(index,None)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for analyze_samples_pca fails'%outfile)
    out_node = bie3.Data(rulebase.PcaAnalysis,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'pca_matrix_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters
    

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

