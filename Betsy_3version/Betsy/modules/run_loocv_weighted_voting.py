#run_loocv.py
import arrayio
import os
import svmutil
from Betsy import bie3,rule_engine_bie3
from Betsy import rulebase
from Betsy import read_label_file
from Betsy import module_utils

def run(in_nodes,parameters, user_input,network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    M = arrayio.read(data_node.identifier)
    a,training_label,second_line = read_label_file.read(
        cls_node.identifier)
    full_index = range(M.ncol())
    predict_model = __import__('Betsy.modules.' + 'classify_with_weighted_voting',globals(),
                                 locals(),['classify_with_weighted_voting'],-2)
    evaluate_model =  __import__('Betsy.modules.' + 'evaluate_prediction',globals(),
                                 locals(),['evaluate_prediction'],-2)
    f = file(outfile,'w')
    f.write('\t'.join(['sample_name','Predicted_class','Confidence','Actual_class','Correct?']))
    f.write('\n')
    for i in range(M.ncol()):
        test_index = i
        train_index = full_index[:]
        train_index.remove(test_index)
        y_training = [training_label[x] for x in train_index]
        y_test = [training_label[test_index]]
        M_train = M.matrix(None,train_index)
        M_test = M.matrix(None,[test_index])
        train_file = 'train'+'_'+str(i)
        test_file = 'test'+'_'+str(i)
        f_out = file(train_file,'w')
        arrayio.gct_format.write(M_train,f_out)
        f_out.close()
        f_out = file(test_file,'w')
        arrayio.gct_format.write(M_test,f_out)
        f_out.close()
        train_label = 'train_label'+'_'+ str(i)
        test_label = 'test_label' + '_' + str(i)
        read_label_file.write(train_label,second_line,y_training)
        read_label_file.write(test_label,second_line,y_test[0])
        train_node = rulebase.PrettySignalFile.output(format='gct',contents='class0,class1')
        train_data = rule_engine_bie3.DataObject(train_node,identifier='train'+'_'+str(i))
        test_node = rulebase.PrettySignalFile.output(format='gct',contents='test')
        test_data = rule_engine_bie3.DataObject(test_node,identifier='test'+'_'+str(i))
        train_label_node = rulebase.ClassLabelFile.output(contents='class0,class1')
        train_label_data = rule_engine_bie3.DataObject(train_label_node,identifier='train_label'+'_'+str(i))
        test_label_node = rulebase.ClassLabelFile.output(contents='test')
        test_label_data = rule_engine_bie3.DataObject(test_label_node,identifier = 'test_label'+'_'+str(i))
        x= train_data,test_data,train_label_data
        out_node = predict_model.run(x,parameters,user_input,network)
        out_node_label=evaluate_model.run((out_node,test_label_data),parameters,user_input,network)
        f1 = open(out_node_label.identifier,'r')
        lines = f1.readlines()
        f1.close()
        f.write(lines[1])
        os.remove(train_file)
        os.remove(test_file)
        os.remove(train_label)
        os.remove(test_label)
        os.remove(out_node.identifier)
        os.remove(out_node_label.identifier)
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
    filename = 'predication_loocv_wv' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)
