#filter_genes_by_fold_change_across_classes.py
import os
#from Betsy
import module_utils, read_label_file
from genomicode import jmath
import math
import bie
import rulebase

def run(x, parameters):
    import arrayio
    data_node,cls_node = x
    outfile = name_outfile(x)
    # obtain the class label
    label, label_line, second_line = read_label_file.read(
        cls_node.attributes['filename'])
    class_num = len(label)
    assert class_num == 2, 'the number of class is not 2'
    fc = int(parameters['group_fc'])
    M = arrayio.read(data_node.attributes['filename'])
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        fold_change = abs(jmath.mean(first[i])-jmath.mean(second[i]))
        if fold_change >= math.log(fc,2):
            I_good.append(i)
    assert I_good, 'there is no gene is significant in fold change with 2'
    f = file(outfile,'w')
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for filter_genes_by_fold_change_across_classes fails'
         % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
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
    filename = 'signal_group_fc_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
