#apply_svm_model.py
import svmutil
import sys
import arrayio
import read_label_file
import rule_engine
import module_utils
import os
import read_label_file
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    model_file,obj = module_utils.find_object(
        parameters,objects,'svm_model','DatasetId')
    assert os.path.exists(model_file),'the model file does not exist'
    test=arrayio.read(identifier)
    x_test = module_utils.format_convert(test)#convert to the format libsvm accept
    test_label_file,obj = module_utils.find_object(parameters,objects,
                                'class_label_file','TestContents,DatasetId')
    assert os.path.exists(test_label_file),'the input file does not exist'
    a,test_label,second_line = read_label_file.read(test_label_file)    
    y_test = [int(x) for x in test_label]
    model = svmutil.svm_load_model(model_file)
    p_label,p_acc,p_val = svmutil.svm_predict(y_test,x_test,model)
    train_label_file,obj = module_utils.find_object(
        parameters,objects,'class_label_file','TrainContents,DatasetId')
    assert os.path.exists(train_label_file),'cannot find train_label_file'
    result,label_line,second_line = read_label_file.read(train_label_file)
    legend_name = ['"'+ i +'='+
                    str(second_line.index(i))+'"' for i in second_line]
    from genomicode import jmath
    R = jmath.start_R()
    R('library(R.utils)')
    command = 'png2("'+outfile+'")'
    R(command)
    jmath.R_equals_vector(p_label,'p_label')
    R('barx <- barplot(p_label, ylim=c(-1,max(p_label)+1), \
      col="blue", axis.lty=2, ylab="Prediction",xlab="Sample")')
    jmath.R_equals_vector(legend_name,'legend_name')
    R('legend("bottomleft", legend_name, lty=1,cex=1)')
    R('dev.off()')
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(parameters,
            objects,'signal_file','TestContents,DatasetId',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','TestContents,DatasetId')
    assert os.path.exists(identifier),'the test file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'classification_file',parameters,objects,single_object)
    return new_objects
    
