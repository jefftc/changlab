#rank_genes_by_sample_ttest.py
#from Betsy
import gene_ranking
from Betsy import module_utils
import shutil
import os
from genomicode import jmath
import arrayio
import numpy
from Betsy import read_label_file,bie,rulebase


def run(in_nodes, parameters, network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    label,label_line,second_line = read_label_file.read(
        cls_node.attributes['filename'])
    M = arrayio.read(data_node.attributes['filename'])
    assert len(label) == 2, (
        'the length of label in %s should be 2'%cls_node.attributes['filename'])
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None,label[0][0])
    second = M.slice(None,label[1][0])
    t,p = gene_ranking.t_test(first,second)
    for i in range(len(p)):
        if not p[i]:
            p[i]=10
    sort_p = [(p[index],index) for index in range(len(p))]
    key = M._row_order[0]
    sort_p.sort()
    gene_list=[]
    key = M._row_order[0]
    threshold = 0.05
    if parameters['gene_select_threshold']:
        threshold = float(parameters['gene_select_threshold'])
    if parameters['gene_order'] == 't_test_p':
        for i in range(len(sort_p)):
            if float(sort_p[i][0]) < threshold:
                gene_list.append(M._row_names[key][sort_p[i][1]])
    elif parameters['gene_order'] == 't_test_fdr':
        for i in range(len(p)):
            if p[i] == 10:
                p[i]= ''
        fdr = jmath.cmh_fdr_bh(p)
        for i in range(len(fdr)):
            if numpy.isnan(fdr[i]):
                fdr[i]=10
        sort_fdr = [(fdr[index],index) for index in range(len(fdr))]
        sort_fdr.sort()
        for i in range(len(fdr)):
            if float(sort_fdr[i][0]) < threshold:
                gene_list.append(M._row_names[key][sort_fdr[i][1]])
    f = open(outfile,'w')
    f.write('\t'.join(gene_list))
    f.close()
    assert len(gene_list)>0,'there is no significant genes can be found in ttest'
    assert module_utils.exists_nz(outfile),(
        'the output file %s for rank_genes_by_sample_ttest fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.GeneListFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'gene_list' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


   
