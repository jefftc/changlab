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
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    training = arrayio.read(identifier)
    x_training = module_utils.format_convert(training)#convert to the format libsvm accept
    training_label_file,obj=module_utils.find_object(parameters,
                                    objects,'class_label_file','Contents,DatsetId')
    a,training_label,second_line = read_label_file.read(training_label_file)
    y_training = [int(x) for x in training_label]
    model = svmutil.svm_train(y_training,x_training)
    svmutil.svm_save_model(outfile,model)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatsetId',pipeline)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatsetId','svm_model',pipeline)

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatsetId')
    



