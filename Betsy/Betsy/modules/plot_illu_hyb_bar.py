#plot_illu_hyb_bar.py
import os
from Betsy import module_utils
import shutil
from genomicode import mplgraph
import numpy
import math
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    plot_hyb_bar(single_object.identifier,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_hyb_bar fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'hyb_bar_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'control_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_ilu_hyb_bar does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'plot_illu_hyb_bar',parameters,objects,single_object)
    return new_objects

def plot_hyb_bar(filename,outfile):
    high = ['ILMN_2038770','ILMN_2038769']
    med = ['ILMN_2038768','ILMN_2038771']
    low = ['ILMN_1343050','ILMN_1343052']
    high_data = []
    med_data = []
    low_data = []
    import arrayio
    M = arrayio.read(filename)
    header = M.row_names()
    for i in range(M.dim()[0]):
        if not M.row_names(header[1])[i] == 'cy3_hyb':
            continue
        if M.row_names(header[0])[i] in high:
            high_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in med:
            med_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in low:
            low_data.extend(M.slice()[i])
    mean=[numpy.mean(high_data),numpy.mean(med_data),numpy.mean(low_data)]
    flag = [math.isnan(i) for i in mean]
    assert True not in flag,'input is not a control file'
    std=[numpy.std(high_data),numpy.std(med_data),numpy.std(low_data)]
    fig=mplgraph.barplot(mean,std,ylabel='Signal',box_label=['high','med','low'])
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile),'the plot_illu_hyb_bar.py fails'
        
