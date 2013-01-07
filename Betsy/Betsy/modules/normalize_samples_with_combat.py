#normalize_samples_with_combat.py
import os
import subprocess
from Betsy import module_utils, read_label_file
from genomicode import config

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
    parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),'cannot find label_file %s'%label_file.identifier
    result,label_line,second_line=read_label_file.read(label_file.identifier)
    assert len(result) >= 2, 'for combat,there should be equal or larger than 2 classes'
    combat_path = config.combatnorm
    combat_BIN = module_utils.which(combat_path)
    assert combat_BIN,'cannot find the %s' %combat_path
    command = ['python', combat_BIN,'-f',single_object.identifier,'-o',
               outfile,'-label',label_file.identifier]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for combat fails' %outfile)
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
    filename = 'signal_combat_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for combat does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
    


