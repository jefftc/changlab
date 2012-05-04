#plot_prediction.py
import sys
import rule_engine
import module_utils
import os
from genomicode import jmath

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    f = open(identifier,'r')
    text = f.readlines()
    f.close()
    label_dict = dict()
    label_list = []
    i = -1
    for line in text[1:]:
        label = line.split()[1]
        if  label not in lable_dict.keys():
            i = i+1
            label_dict[label] = i
        label_list.append(label_dict[label])
    R = jmath.start_R()
    R('library(R.utils)')
    command = 'png2("'+outfile+'")'
    R(command)
    legend_name = ['"'+label_dict[i] + '=' +
                    i+'"' for i in label_dict.keys()]
    jmath.R_equals_vector(label_list,'p_label')
    jmath.R_equals_vector(legend_name,'legend_name')
    R('barx <- barplot(p_label, ylim=c(-1,max(p_label)+1), \
      col="blue", axis.lty=2, ylab="Prediction",xlab="Sample")')
    R('legend("bottomleft", legend_name, lty=1,cex=1)')
    R('dev.off()')
    assert module_utils.exists_nz(outfile),'the output\
                    file %s for plot_prediction fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(parameters,
            objects,'prediction_plot','testcontents',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'classification_file','traincontents')
    assert os.path.exists(identifier),'the classification\
                file %s for plot_prediction does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'classification_file',parameters,objects,single_object)
    return new_objects
    



