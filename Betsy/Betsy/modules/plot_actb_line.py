#plot_actb_line.py

import os
from Betsy import module_utils
import shutil
from genomicode import mplgraph,arrayplatformlib, config
import arrayio
import subprocess
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    header = M._row_names
    header_name = M._row_order
    column_id = None
    for i in header_name:
        column = [j.upper() for j in header[i]]
        if 'ACTB' in column:
            column_id = i
            break
    if not column_id:
        x = arrayplatformlib.identify_all_platforms_of_matrix(M)
        id = x[0][0]
        platform = x[0][1]
        if platform in ['HG_U133_Plus_2','HG_U133B','HG_U133A','HG_U133A_2','HG_U95A',
                        'HumanHT_12','HG_U95Av2','entrez_ID_human','entrez_symbol_human','Hu6800']:
            out_platform = 'Entrez_symbol_human'
        elif platform in ['Mouse430A_2','MG_U74Cv2', 'Mu11KsubB','Mu11KsubA','MG_U74Av2',
             'Mouse430_2', 'MG_U74Bv2','entrez_ID_mouse','MouseRef_8','entrez_symbol_mouse']:
            out_platform= 'Entrez_symbol_mouse'
        Annot_path = config.annotate_matrix
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN,'cannot find the %s' %Annot_path
        command = ['python', Annot_BIN,'-f',single_object.identifier,'-o','tmp','--platform',
               out_platform]
        process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        assert module_utils.exists_nz('tmp'),'the platform conversion fails'
        column_id = out_platform
        M = arrayio.read('tmp')
    label = M._col_names['_SAMPLE_NAME']
    keywords = ['ACTB','TUBB']
    lines= []
    data=[]
    legend_name = []
    for keyword in keywords:
        for i in range(M.dim()[0]):
            if M._row_names[column_id][i].upper() == keyword:
                data.append(M.slice()[i])
                legend_name.append(keyword+'('+M._row_names[header_name[0]][i]+')')
    for i in range(len(data)):
        line = [(j,data[i][j]) for j in range(len(data[i]))]
        lines.append(line)
        fig=mplgraph.lineplot(*lines,legend=legend_name,box_label=label,
                          ylim_min=0,ylabel='Gene Expression Value')
        fig.savefig(outfile)
    if not lines:
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot([0,0,0,0])
        plt.title('no ACTB or TUBB probes are found')
        plt.savefig(outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_actb_line fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                      parameters,single_object,pipeline,outfile,starttime,user,jobname)
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
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_actb_line does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'actb_plot',parameters,objects,single_object)
    return new_objects
