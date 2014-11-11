#classify_with_random_forest.py

import sys
import arrayio
import os
from Betsy import read_label_file, module_utils
from genomicode import jmath
from Betsy import bie3, rulebase

def run(in_nodes, parameters, user_input, network,num_cores):
    data_node,cls_node_train = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    result, label_line, second_line = read_label_file.read(
       cls_node_train.identifier)
    y=[second_line[int(i)] for i in label_line]
    R = jmath.start_R()
    M = arrayio.read(data_node.identifier)
    M_train = M.matrix(None,range(0,len(label_line)))
    M_test = M.matrix(None,range(len(label_line),M.dim()[1]))
    M1 = M_train.slice()
    M_train = jmath.transpose(M1)
    jmath.R_equals_matrix(M_train, 'data')
    M2 = M_test.slice()
    M2 = jmath.transpose(M2)
    jmath.R_equals_matrix(M2, 'test')
    jmath.R_equals(y, 'y')
    R('y<-as.factor(y)')
    R('require(randomForest,quietly=TRUE)')
    R('library(randomForest)')
    R('model<-randomForest(data,y=y,importance=TRUE)')
    R('predict_result<-predict(model, test)')
    predict_result = R['predict_result']
    levels = predict_result.levels
    predict_labels = predict_result[:]
    predict_labels = [levels[i - 1] for i in predict_labels]
    name = M_test._col_names.keys()[0]
    sample_name = M_test._col_names[name]
    result = [['Sample_name', 'Predicted_class', 'Confidence']]
    for i in range(len(sample_name)):
        result.append(
            [str(sample_name[i]), predict_labels[i], ''])
    f = file(outfile, 'w')
    for i in result:
        f.write('\t'.join(i))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for classify_with_random_forest fails' % outfile)
    out_node = bie3.Data(rulebase.ClassifyFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def name_outfile(in_nodes,user_input):
    data_node,cls_node_train = in_nodes
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'random_forest_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node_train = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            contents='class0,class1,test',
                                            datatype='SignalFile')
    
    cls_node_train = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                                 contents='class0,class1',
                                            datatype='ClassLabelFile')
    
    return data_node,cls_node_train
