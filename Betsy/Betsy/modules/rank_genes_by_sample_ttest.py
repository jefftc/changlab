#rank_genes_by_sample_ttest.py
from Betsy import gene_ranking
from Betsy import module_utils
import os
from genomicode import jmath
import arrayio
import numpy
from Betsy import read_label_file, bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    label, label_line, second_line = read_label_file.read(cls_node.identifier)
    M = arrayio.read(data_node.identifier)
    assert len(label) == 2, (
        'the length of label in %s should be 2' % cls_node.identifier
    )
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    t, p = gene_ranking.t_test(first, second)
    for i in range(len(p)):
        if not p[i]:
            p[i] = 10
    sort_p = [(p[index], index) for index in range(len(p))]
    key = M._row_order[0]
    sort_p.sort()
    gene_list = []
    key = M._row_order[0]
    threshold = 0.05
    if 'gene_select_threshold' in user_options:
        threshold = float(user_options['gene_select_threshold'])
    if out_attributes['gene_order'] == 't_test_p':
        for i in range(len(sort_p)):
            if float(sort_p[i][0]) < threshold:
                gene_list.append(M._row_names[key][sort_p[i][1]])
    elif out_attributes['gene_order'] == 't_test_fdr':
        for i in range(len(p)):
            if p[i] == 10:
                p[i] = ''
        fdr = jmath.cmh_fdr_bh(p)
        for i in range(len(fdr)):
            if numpy.isnan(fdr[i]):
                fdr[i] = 10
        sort_fdr = [(fdr[index], index) for index in range(len(fdr))]
        sort_fdr.sort()
        for i in range(len(fdr)):
            if float(sort_fdr[i][0]) < threshold:
                gene_list.append(M._row_names[key][sort_fdr[i][1]])
    f = open(outfile, 'w')
    f.write('\t'.join(gene_list))
    f.close()
    assert len(gene_list) > 0, 'there is no significant genes can be found in ttest'
    assert module_utils.exists_nz(outfile), (
        'the output file %s for rank_genes_by_sample_ttest fails' % outfile
    )
    out_node = bie3.Data(rulebase.GeneListFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='_SignalFile_Order')
    cls_node = module_utils.get_identifier(network, module_id, pool,
                                           user_attributes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'gene_list' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
