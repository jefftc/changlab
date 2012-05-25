#run_scoresig.py
import os
import hash_method
import module_utils
import rule_engine
import subprocess
def run(parameters,objects,pipeline):
    rma_file,obj1 = module_utils.find_object(parameters,objects,'signal_file','contents,pre1')
    mas_file,obj2 = module_utils.find_object(parameters,objects,'signal_file','contents,pre2')
    assert os.path.exists(rma_file),'the rma_file %s in run_scoresig does not exist'%rma_file
    assert os.path.exists(mas_file),'the mas_file %s in run_scoresig does not exist'%mas_file
    outfile = get_outfile(parameters,objects,pipeline)
    import Betsy_config
    scoresig_path = Betsy_config.SCORESIG
    scoresig_BIN = module_utils.which(scoresig_path)
    assert scoresig_BIN,'cannot find the %s' %scoresig_path

    command = ['python',scoresig_BIN,'-r',rma_file,'-m',mas_file,'-j','20','-o',outfile]
    process = subprocess.Popen(command,shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),'the output \
        file %s for run_scoresig does not exists'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                            obj1,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)


def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
    parameters,objects,'signal_file','contents',pipeline)
    

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input file %s for run_scoresig does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signature_score',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
