#plot_signal_heatmap.py
import os
import subprocess
from Betsy import module_utils
from genomicode import config, graphlib
import arrayio
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    """generate a heatmap of input file"""
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    Heatmap_path = config.arrayplot
    Heatmap_BIN = module_utils.which(Heatmap_path)
    assert Heatmap_BIN,'cannot find the %s' %Heatmap_path

    command = ['python', Heatmap_BIN,single_object.identifier,'-o',outfile,"--label_arrays",
               "--label_genes",'--no_autoscale']#"--grid"
    if 'color' in parameters.keys():
        color=['--color' , parameters['color'].replace('_','-')]
        command.extend(color)
    M = arrayio.read(single_object.identifier)
    nrow = M.nrow()
    ncol = M.ncol()
    ratio = float(nrow)/ncol
    max_box_height = None
    max_box_width = None
    if 'hm_width' in parameters.keys():
        max_box_width = parameters['hm_width']
    if 'hm_height' in parameters.keys():
         max_box_height = parameters['hm_height']
    if ratio >= 4:
        x,y=graphlib.find_tall_heatmap_size(nrow,ncol,
                                            max_box_height=max_box_height,
                                            max_box_width=max_box_width,
                                            max_megapixels=128)
    else:
        x,y=graphlib.find_wide_heatmap_size(nrow,ncol,
                                            max_box_height=max_box_height,
                                            max_box_width=max_box_width,
                                            max_megapixels=128)
    command.extend(['-x',str(x),'-y',str(y)])
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_signal_heatmap fails' %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile,starttime,user,jobname)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
            identifier,pipeline,parameters)
   
def get_identifier(parameters,objects):
    if parameters['cluster_alg'] == 'no_cluster_alg':
        single_object = module_utils.find_object(
            parameters,objects,'signal_file','contents')
    else:
        single_object = module_utils.find_object(
            parameters,objects,'cluster_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_signal_heatmap does not exist'
        %single_object.identifier)
    return single_object

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'heatmap'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'cluster_heatmap',parameters,objects,single_object)
    return new_objects
    
