#filter_and_threshold_genes.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase

def run(data_node, parameters, user_input, network):
    """run preprocessdataset """
    import arrayio
    outfile = name_outfile(data_node,user_input)
    threshold = 20
    ceiling = 16000
    min_fold_change = 5
    min_delta = 100.0
    M = arrayio.read(data_node.identifier)
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        for j in range(len(X[i])):
            if X[i][j]<threshold:
                M._X[i][j] = threshold
            if X[i][j]>ceiling:
                M._X[i][j]=ceiling
        gene = M._X[i]
        fold_change = max(gene)/float(min(gene))
        delta = max(gene)-min(gene)
        if fold_change >= min_fold_change and delta>=min_delta:
            I_good.append(i)
    f = file(outfile,'w')
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for filter_and_threshold_genes fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_preprocessdataset_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


