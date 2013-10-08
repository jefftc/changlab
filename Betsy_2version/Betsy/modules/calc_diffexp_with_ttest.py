#calc_diffexp_with_ttest.py
import shutil
import os
import arrayio
import numpy
from genomicode import jmath
#from Betsy
import gene_ranking, module_utils, read_label_file
from time import strftime,localtime
import bie
import rulebase

def run(in_nodes, parameters):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    label, label_line, second_line = read_label_file.read(
        cls_node.attributes['filename'])
    M = arrayio.read(data_node.attributes['filename'])
    assert len(label) == 2, (
        'the length of label in %s should be 2' % cls_node.attributes['filename'])
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    t, p = gene_ranking.t_test(first, second)
    higher_group = get_higherexpression(M, label, second_line)
    bonf = jmath.cmh_bonferroni(p)
    fdr = jmath.cmh_fdr_bh(p)
    p_copy = p[:]
    for i in range(len(p_copy)):
        if not p_copy[i]:
            p_copy[i] = 10
    sort_p = [(p_copy[index], index) for index in range(len(p_copy))]
    sort_p.sort()
    f = file(outfile, 'w')
    header = ['#']
    header.extend(M._row_order[:])
    header.extend(['p_value', 'cmh_bonferroni',
                   'cmh_fdr', 'higher_expression'])
    f.write('\t'.join(header))
    f.write('\n')
    p = [' ' if not x else x for x in p]
    bonf = [' ' if not x else x for x in bonf]
    fdr = [' ' if not x else x for x in fdr]
    for i in range(len(p_copy)):
        f.write(str(i + 1) + '\t')
        for key in M._row_order:
            f.write(M._row_names[key][sort_p[i][1]])
            f.write('\t')
        f.write(str(p[sort_p[i][1]]) + '\t')
        f.write(str(bonf[sort_p[i][1]]) + '\t')
        f.write(str(fdr[sort_p[i][1]]) + '\t')
        f.write(str(higher_group[sort_p[i][1]]) + '\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for t_test fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.DiffExprFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 't_test_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def get_higherexpression(M, label, second_line):
    higher_group = []
    assert len(label) == 2, 'the length of label should be 2'
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    for i in range(M.nrow()):
        group1 = sum(first[i]) / float(len(first[i]))
        group2 = sum(second[i]) / float(len(second[i]))
        if group1 >= group2:
            higher_group.append(second_line[int(label[0][1])])
        else:
            higher_group.append(second_line[int(label[1][1])])
    return higher_group
