#plot_actb.py

import os
import module_utils
import shutil
from genomicode import mplgraph
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    header = M.row_names()
    label = M._col_names['_SAMPLE_NAME']
    keywords = ['ACTB','TUBB']
    lines= []
    data=[]
    legend_name = []
    for keyword in keywords:
        for i in range(M.dim()[0]):
            if M.row_names(header[1])[i] == keyword:
                data.append(M.slice()[i])
                legend_name.append(keyword+'('+M.row_names(header[0])[i]+')')
        assert len(data)>0,'input is not a control file'
    for i in range(len(data)):
        line = [(j,data[i][j]) for j in range(len(data[i]))]
        lines.append(line)
    fig=mplgraph.lineplot(*lines,legend=legend_name,box_label=label,
                          ylim_min=0,ylabel='Gene Expression Value')
    fig.savefig(outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_actb fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                      parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'actb_plot_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_actb does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'actb_plot',parameters,objects,single_object)
    return new_objects
