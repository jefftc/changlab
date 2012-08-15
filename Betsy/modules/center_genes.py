#center_genes.py
import os
import subprocess
import module_utils

def run(parameters,objects,pipeline):
    """mean or median"""
    CLUSTER_BIN = 'cluster'
    center_alg = {'mean':'a','median':'m'}
    try :
        center_parameter = center_alg[parameters['gene_center']]
    except:
        raise ValueError("Centering parameter is not recognized")
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    
    process = subprocess.Popen([CLUSTER_BIN,'-f',single_object.identifier,
                                '-cg',center_parameter,'-u',outfile],
                                shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    outputfile = outfile + '.nrm'
    os.rename(outputfile,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for centering fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_centering_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for centering does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
