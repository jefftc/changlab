#analyze_samples_pca.py
import os
from genomicode import pcalib
import arrayio
#from Betsy
import module_utils
from time import strftime,localtime
import bie
import rulebase

def run(data_node,parameters):
    outfile = name_outfile(data_node)
    M = arrayio.read(data_node.attributes['filename'])
    X = M._X
    N = int(parameters['pca_gene_num'])
    N = min(N,M.nrow())
    index = pcalib.select_genes_var(X,N)
    M_new = M.matrix(index,None)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for analyze_samples_pca fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.PcaAnalysis,**new_parameters)
    return out_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'pca_matrix_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters
    

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

