#plot_sample_pca_with_predictions.py
import os
from Betsy import module_utils
import shutil
from Betsy import read_label_file
from genomicode import pcalib,genesetlib
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    legend_object = module_utils.find_object(
        parameters,objects,'classification_file','traincontents')
    assert os.path.exists(legend_object.identifier),(
        'the prediction file %s for plot_sample_pca_with_predictions does not exist'
        %legend_object.identifier)
    result_data = genesetlib.read_tdf(legend_object.identifier,
                                      preserve_spaces=True,allow_duplicates=True)
    for i in result_data:
        if i[0] == 'Predicted_class':
            legend = i[2]
    colors = ['r','b','g','y']
    legend_dict = {}
    for index,item in enumerate(legend):
        if item not in legend_dict:
            legend_dict[item]=[index]
        else:
            legend_dict[item].append(index)
    color=['']*len(legend)
    for index,key in enumerate(legend_dict.keys()):
        c = colors[index]
        for i in legend_dict[key]:
            color[i]=c
    module_utils.plot_pca(single_object.identifier,outfile,color,legend)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_sample_pca_with_predictions fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    filename = 'prediction_pca_plot.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'pca_matrix','testcontents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_sample_pca_with_predictions does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'prediction_pca_plot',parameters,objects,single_object)
    return new_objects



