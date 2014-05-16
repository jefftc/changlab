#filter_genes_by_missing_values.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase

def run(data_node,parameters, user_input, network):
    outfile = name_outfile(data_node,user_input)
    import arrayio
    f_out = file(outfile,'w')
    M = arrayio.read(data_node.identifier)
    I_good = []
    #get the percentage of gene filter
    percent=float(user_input['filter_value'])/100
    for i in range(M.dim()[0]):
       missing_count = 0
       for j in range(M.dim()[1]):
           if M._X[i][j] in [None,'NA']:
                missing_count = missing_count + 1
       if float(missing_count)/M.dim()[1]<percent:
            I_good.append(i)
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for gene_filter fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile_Impute,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_filter_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

