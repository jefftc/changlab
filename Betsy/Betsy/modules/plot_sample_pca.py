#plot_sample_pca.py
import os
from Betsy import module_utils
import shutil
from Betsy import read_label_file
import arrayio
import matplotlib.cm as cm
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    a,b,c = read_label_file.read(label_file.identifier)
    if len(a)>1:
        colors = []
        for i in range(5):
            colors.append(cm.hot(i/5.0,1))
            colors.append(cm.autumn(i/5.0,i))
            colors.append(cm.cool(i/5.0,i))
            colors.append(cm.jet(i/5.0,i))
            colors.append(cm.spring(i/5.0,i))
            colors.append(cm.prism(i/5.0,i))
            colors.append(cm.summer(i/5.0,i))
            colors.append(cm.winter(i/5.0,i))
        opts = [colors[int(i)] for i in b]
        legend = [c[int(i)] for i in b]
        module_utils.plot_pca(single_object.identifier,outfile,opts,legend)
    else:
        module_utils.plot_pca(single_object.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for pca_sample_plot fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                             [single_object,label_file],
                                             pipeline,outfile,starttime,user,jobname)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = parameters['objecttype']+'_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'pca_matrix','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for pca_sample_plot does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'pca_plot',parameters,objects,single_object)
    return new_objects



