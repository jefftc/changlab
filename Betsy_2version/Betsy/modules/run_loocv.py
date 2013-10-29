#run_loocv.py
import arrayio
from Betsy import read_label_file
from Betsy import module_utils, bie, rulebase
import os
import svmutil
 

def run(in_nodes,parameters, network):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    M = arrayio.read(data_node.attributes['filename'])
    a,training_label,second_line = read_label_file.read(cls_node.attributes['filename'])
    full_index = range(M.ncol())
    train_model = None
    model_dict = {'svm':'classify_with_svm','svm_test':'train_svm_model',
                  'random_forest':'classify_with_random_forest',
                  'weighted_voting':'classify_with_weighted_voting'}
    train_model = None
    if parameters['classify_alg'] == 'svm':
        train_model = __import__('modules.' + model_dict['svm_test'],globals(),
                                 locals(),[model_dict['svm_test']],-2)
    predict_model = __import__('modules.' + model_dict[parameters['classify_alg']],globals(),
                                 locals(),[model_dict[parameters['classify_alg']]],-2)
        
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
        train_node = SignalFile2(format='gct',contents='class0,class1',filename='train'+'_'+str(i))
        test_node = SignalFile2(format='gct',contents='test',filename='test'+'_'+str(i))
        train_label_node = ClassLabelFile(contents='class0,class1',filename='train_label'+'_'+str(i))
        test_label_node = ClassLabelFile(contents='test',filename='test_label'+'_'+str(i))
        if train_model:
            new_parameters = {}
            for key in SvmModel.get_defaults():
                new_parameters[key]=parameters[key]
            x1 = train_node,train_label_node
            svm_model = train_model.run(x1,new_parameters)
            x = svm_model,test_node,test_label_node
        else:
            x= train_node,test_node,train_label_node,test_label_node
        out_node = predict_model.run(x,parameters)
        f1 = open(out_node.attributes['filename'],'r')
        lines = f1.readlines()
        f1.close()
        f.write(lines[1])
        os.remove(train_file)
        os.remove(test_file)
        os.remove(train_label)
        os.remove(test_label)
        os.remove(out_node.attributes['filename'])
        if train_model:
            os.remove(svm_model.attributes['filename'])
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for loocv fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.ClassifyFile,**new_parameters)
    return out_node
    


def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2',
                                            contents='class0,class1')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile',
                                           contents='class0,class1')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'predication_loocv_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
