#plot_intensity.py
import os
import module_utils
import shutil
import math
from genomicode import mplgraph,jmath
import arrayio
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    data = jmath.transpose(M._X)
    tickname=M._col_names['_SAMPLE_NAME']
    fig = mplgraph.boxplot(data,xlabel='sample_name',ylabel='signal',title='signal_intensity',box_label=tickname)
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot intensity fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'intensity_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_intensity does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'intensity_plot',parameters,objects,single_object)
    return new_objects

