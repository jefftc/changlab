#train_svm_model.py
import svmutil
import sys
import arrayio
from Betsy import read_label_file
from Betsy import module_utils,bie,rulebase
import os


def run(in_nodes, parameters, network):
    data_node_train,cls_node_train = in_nodes
    outfile = name_outfile(in_nodes)
    a,training_label,second_line = read_label_file.read(
        cls_node_train.attributes['filename'])
    training = arrayio.read(data_node_train.attributes['filename'])
    x_training = module_utils.format_convert(training.matrix(None,range(0,len(training_label))))#convert to the format libsvm accept
    y_training = [int(x) for x in training_label]
    svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel']
    #if 'svm_kernel' in parameters.keys():
    kernel_type = svm_kernel.index(parameters['svm_kernel'])
    command = '-t '+ str(kernel_type)
   # else:
    #    command = '-t 0'
    param = svmutil.svm_parameter(command)
    prob  = svmutil.svm_problem(y_training, x_training)
    model = svmutil.svm_train(prob,param)
    svmutil.svm_save_model(outfile,model)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for train_svm_model fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SvmModel,**new_parameters)
    return out_node
    

def name_outfile(in_nodes):
    data_node_train,cls_node_train= in_nodes
    original_file = module_utils.get_inputid(data_node_train.attributes['filename'])
    filename = 'svm_model_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node_train,cls_node_train = in_nodes
    identifier = data_node_train.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node_train = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='class0,class1,test',
                                            datatype='SignalFile2')
    cls_node_train = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='class0,class1',
                                            datatype='ClassLabelFile')
    return data_node_train,cls_node_train
