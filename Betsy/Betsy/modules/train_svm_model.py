#train_svm_model.py
import svmutil
import sys
import arrayio
from Betsy import read_label_file
from Betsy import module_utils
import os
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    training = arrayio.read(single_object.identifier)
    x_training = module_utils.format_convert(training)#convert to the format libsvm accept
    training_label_file = module_utils.find_object(parameters,
                                    objects,'class_label_file','contents')
    assert os.path.exists(training_label_file.identifier),(
        'the training label file %s does not exist'
        %training_label_file.identifier)
    a,training_label,second_line = read_label_file.read(training_label_file.identifier)
    y_training = [int(x) for x in training_label]
    svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel']
    if 'svm_kernel' in parameters.keys():
        kernel_type = svm_kernel.index(parameters['svm_kernel'])
        command = '-t '+ str(kernel_type)
    else:
        command = '-t 0'
    param = svmutil.svm_parameter(command)
    prob  = svmutil.svm_problem(y_training, x_training)
    model = svmutil.svm_train(prob,param)
    svmutil.svm_save_model(outfile,model)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for train_svm_model fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'svm_model_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
                           parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for train_svm_model does not exist'
        %single_object.identifier)
    return single_object


def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = module_utils.DataObject('svm_model',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
