#generate_genelist_from_diffexprfile.py
from Betsy import module_utils
import os
import math
from Betsy import read_label_file,bie3,rulebase


def run(data_node, parameters, user_input, network):
    outfile = name_outfile(data_node,user_input)
    foldchange = None
    p_value = None
    fdr = None
    foldchange_index = None
    p_value_index = None
    fdr_index = None
    if 'select_gene_by_foldchange' in user_input and user_input['select_gene_by_foldchange']:
        foldchange = float(user_input['select_gene_by_foldchange'])
    if 'select_gene_by_p_value' in user_input and user_input['select_gene_by_p_value']:
        p_value = float(user_input['select_gene_by_p_value'])
    if 'select_gene_by_fdr' in user_input and user_input['select_gene_by_fdr']:
        fdr = float(user_input['select_gene_by_fdr'])
    
    f = file(data_node.identifier,'r')
    text = f.readlines()
    f.close()
    gene_list = range(1,len(text))
    headers = text[0].split('\t')
    if foldchange:
        new_gene_list = []
        foldchange_index = headers.index('Log_2 Fold Change')
        for index in gene_list:
            line = text[index]
            values = line.split('\t')
            fc = values[foldchange_index]
            if fc>=math.log(foldchange,2):
                new_gene_list.append(index)
        gene_list = new_gene_list
    if p_value:
        p_value_index = headers.index('NL10P')
        new_gene_list = []
        for index in gene_list:
            line = text[index]
            values = line.split('\t')
            pv = values[p_value_index]
            if pv>=math.log(p_value,10):
                new_gene_list.append(index)
        gene_list = new_gene_list
        
    if fdr:
        fdr_index = headers.index('NL10 FDR')
        new_gene_list = []
        for index in gene_list:
            line = text[index]
            values = line.split('\t')
            fdrv = values[fdr_index]
            if fdrv>=math.log(fdr,10):
                new_gene_list.append(index)
        gene_list = new_gene_list
    
    gene_index = headers.index('Gene ID')
    gene_list_name = [text[i].split('\t')[gene_index] for i in gene_list]
            
    f = open(outfile,'w')
    f.write('\t'.join(gene_list_name))
    f.close()
    assert len(gene_list_name)>0,'there is no genes can be found'
    assert module_utils.exists_nz(outfile),(
        'the output file %s for generate_genelist_from_diffexprfile fails'%outfile)
    out_node = bie3.Data(rulebase.GeneListFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='DiffExprFile')
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'gene_list' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    return module_utils.make_unique_hash(data_node.identifier,
                                         pipeline,parameters,user_input)


   
