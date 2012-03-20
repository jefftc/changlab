#train_svm_model.py
import svmutil
import sys
import arrayio
import read_label_file
import rule_engine
import module_utils
import os
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    training = arrayio.read(identifier)
    x_training = module_utils.format_convert(training)#convert to the format libsvm accept
    training_label_file,obj = module_utils.find_object(parameters,
                                    objects,'class_label_file','Contents,DatasetId')
    assert os.path.exists(training_label_file),'the training label file does not exist'
    a,training_label,second_line = read_label_file.read(training_label_file)
    y_training = [int(x) for x in training_label]
    model = svmutil.svm_train(y_training,x_training)
    svmutil.svm_save_model(outfile,model)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object


def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'svm_model',parameters,objects,single_object)
    return new_objects
