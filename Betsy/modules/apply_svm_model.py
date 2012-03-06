#apply_svm_model.py
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
    model_file,obj = module_utils.find_object(
        parameters,objects,'svm_model','Train_DatasetId')
    test=arrayio.read(identifier)
    x_test = module_utils.format_convert(test)#convert to the format libsvm accept
    test_label_file,obj = module_utils.find_object(parameters,objects,'class_label_file','TestContents,Test_DatasetId')
    a,test_label,second_line = read_label_file.read(test_label_file)    
    if test_label:
        y_test = [int(x) for x in test_label]
    else:
        y_test = [0]*len(x_test)
    model = svmutil.svm_load_model(model_file)
    p_label,p_acc,p_val = svmutil.svm_predict(y_test,x_test,model)
    f = file(outfile,'w')
    f.write('\t'.join(str(p_label))+'\n')
    if test_label:
        f.write('\t'.join(str(p_acc)))
        f.write('\t'.join(str(p_val)))
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','TestContents,Test_DatasetId',pipeline)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(parameters,
            objects,'signal_file','TestContents,Test_DatasetId','classification_file',pipeline)

def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'signal_file','TestContents,Test_DatasetId')


