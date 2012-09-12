#cluster_heatmap.py
import os
import hash_method
import subprocess
import module_utils

def run(parameters,objects,pipeline,options=None):
    """generate a heatmap of input file"""
    single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    import Betsy_config
    Heatmap_path = Betsy_config.ARRAYPLOT
    Heatmap_BIN = module_utils.which(Heatmap_path)
    assert Heatmap_BIN,'cannot find the %s' %Heatmap_path
    command = ['python', Heatmap_BIN,single_object.identifier,'-o',outfile,"--label_arrays",
               "--grid","--label_genes"]
    if 'color' in parameters.keys():
        color=['--color' , parameters['color'].replace('_','-')]
        command.extend(color)
    if 'hm_width' in parameters.keys():
        command.extend(['-x',parameters['hm_width']])
    if 'hm_height' in parameters.keys():
        command.extend(['-y',parameters['hm_height']])
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for cluster_heatmap fails' %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
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
        'the input file %s for cluster_heatmap does not exist'
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
    
