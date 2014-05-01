#run_loocv.py
import arrayio
import os
import svmutil
from Betsy import bie3,rule_engine_bie3
from Betsy import rulebase
from Betsy import read_label_file
from Betsy import module_utils

def run(in_nodes,parameters, user_input, network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    M = arrayio.read(data_node.identifier)
    a,training_label,second_line = read_label_file.read(cls_node.identifier)
    full_index = range(M.ncol())
    predict_model = __import__('Betsy.modules.' + 'classify_with_svm',globals(),
                                 locals(),['classify_with_svm'],-2)
    train_model = __import__('Betsy.modules.' + 'train_svm_model',globals(),
                                 locals(),['train_svm_model'],-2)
    evaluate_model =  __import__('Betsy.modules.' + 'evaluate_prediction',globals(),
                                 locals(),['evaluate_prediction'],-2)  
    f = file(outfile,'w')
    f.write('\t'.join(['sample_name','Predicted_class','Confidence','Actual_class','Correct?']))
    f.write('\n')
    for i in range(M.ncol()):
        test_index = i
        train_index = full_index[:]
        train_index.remove(test_index)
        merge_index = train_index + [test_index]
        y_training = [training_label[x] for x in train_index]
        y_test = [training_label[test_index]]
        M_merge = M.matrix(None,merge_index)
        merge_file = 'merge'+'_'+str(i)
        f_out = file(merge_file,'w')
        arrayio.gct_format.write(M_merge,f_out)
        f_out.close()
        train_label = 'train_label'+'_'+ str(i)
        test_label = 'test_label' + '_' + str(i)
        read_label_file.write(train_label,second_line,y_training)
        read_label_file.write(test_label,second_line,y_test[0])
        merge_node =rulebase.PrettySignalFile.output(format='gct',
                                              contents='class0,class1,test')
        merge_data = rule_engine_bie3.DataObject(merge_node,identifier=merge_file)
        train_label_node = rulebase.ClassLabelFile.output(contents='class0,class1')
        train_label_data = rule_engine_bie3.DataObject(train_label_node,identifier=train_label)
        test_label_node = rulebase.ClassLabelFile.output(contents='test')
        test_label_data = rule_engine_bie3.DataObject(test_label_node,identifier=test_label)
        new_parameters = parameters.copy()
        del new_parameters['loocv']
        del new_parameters['actual_label']
        del new_parameters['wv_feature_stat']
        x1 = merge_data,train_label_data
        svm_model = train_model.run(x1,new_parameters,user_input,network)
        x = svm_model,merge_data,train_label_data
        out_node = predict_model.run(x,parameters,user_input,network)
        out_node_label=evaluate_model.run((out_node,test_label_data),parameters,user_input,network)
        f1 = open(out_node_label.identifier,'r')
        lines = f1.readlines()
        f1.close()
        f.write(lines[1])
        os.remove(merge_file)
        os.remove(train_label)
        os.remove(test_label)
        os.remove(out_node.identifier)
        os.remove(out_node_label.identifier)
        os.remove(svm_model.identifier)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for loocv fails'%outfile)
    out_node = bie3.Data(rulebase.ClassifyFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='PrettySignalFile',
                                            contents='class0,class1')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile',
                                           contents='class0,class1')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'predication_loocv_svm' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)
