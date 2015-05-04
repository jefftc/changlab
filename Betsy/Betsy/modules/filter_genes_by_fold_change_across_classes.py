#filter_genes_by_fold_change_across_classes.py
import os
from Betsy import module_utils, read_label_file
from genomicode import jmath
import math
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    import arrayio
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    # obtain the class label
    label, label_line, second_line = read_label_file.read(cls_node.identifier)
    class_num = len(label)
    assert class_num == 2, 'the number of class is not 2'
    fc = 1
    if 'group_fc_num' in user_options:
        fc = int(user_options['group_fc_num'])
    M = arrayio.read(data_node.identifier)
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        fold_change = abs(jmath.mean(first[i]) - jmath.mean(second[i]))
        if fold_change >= math.log(fc, 2):
            I_good.append(i)
    assert I_good, 'there is no gene is significant in fold change with 2'
    f = file(outfile, 'w')
    M_c = M.matrix(I_good, None)
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for filter_genes_by_fold_change_across_classes fails'
        % outfile)
    out_node = bie3.Data(rulebase._SignalFile_Filter, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='_SignalFile_Filter')
    cls_node = module_utils.get_identifier(network, module_id, pool,
                                           user_attributes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_group_fc_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
