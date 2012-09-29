#score_pathway_with_scoresig.py
import os
from Betsy import module_utils
import subprocess
def run(parameters, objects, pipeline):
    rma_file = module_utils.find_object(parameters,objects, 'signal_file', 'contents,pre1')
    mas_file = module_utils.find_object(parameters,objects, 'signal_file', 'contents,pre2')
    assert os.path.exists(rma_file.identifier), (
        'the rma_file %s in run_scoresig does not exist'% rma_file.identifier)
    assert os.path.exists(mas_file.identifier), (
        'the mas_file %s in run_scoresig does not exist'% mas_file.identifier)
    outfile = get_outfile(parameters, objects, pipeline)
    import config
    scoresig_path = config.SCORESIG
    scoresig_BIN = module_utils.which(scoresig_path)
    assert scoresig_BIN,'cannot find the %s' %scoresig_path
    file1,file2 = module_utils.convert_to_same_platform(rma_file.identifier,mas_file.identifier)
    command = ['python', scoresig_BIN, '-r', file1, '-m', file2, '-j', '20', '-o', outfile]
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for run_scoresig does not exists'% outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                            rma_file, pipeline, outfile)
    return new_objects
    
def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signature_score_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile
    

def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for run_scoresig does not exist'
        % single_object.identifier)
    return single_object

def get_newobjects(parameters,objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    attributes = parameters.values()
    new_object = module_utils.DataObject('signature_score', attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
