#plot_prediction_bar.py
import sys
from Betsy import module_utils
import os
from genomicode import mplgraph,filelib

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    matrix=[x for x in filelib.read_cols(filename)]
    header = matrix[0]
    index = header.index('Confidence')
    matrix=matrix[1:]
    confidence = [float(i[index]) for i in matrix]
    sample=[i[0] for i in matrix]
    if confidence==['']*len(matrix):
        index = header.index('Predicted_class')
        class_value = [i[index] for i in matrix]
        label_dict = dict()
        label_list = []
        i = -1
        for label in class_value:
            if  label not in label_dict.keys():
                i = i+1
                label_dict[label] = i
            label_list.append(label_dict[label])
        yticks=label_dict.keys()
        ytick_pos=[label_dict[i] for i in label_dict.keys()]
        fig=mplgraph.barplot(label_list,box_label=sample,
                   ylim=(-0.5,1.5),ytick_pos=ytick_pos,yticks=yticks,
                   xtick_rotation='vertical',ylabel='Prediction',xlabel='Sample')
        fig.savefig(outfile)
    else:
        fig=mplgraph.barplot(confidence,box_label=sample,
                   ylim=(-1.5,1.5),
                   xtick_rotation='vertical',ylabel='Prediction',xlabel='Sample')
        fig.savefig(outfile)
    
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_prediction_bar fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'prediction' + original_file + '.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'classification_file','traincontents')
    assert os.path.exists(single_object.identifier),(
        'the classification file %s for plot_prediction_bar does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'prediction_plot',parameters,objects,single_object)
    return new_objects
    



