#preprocess.py
import os
import module_utils
import subprocess
import hash_method
import rule_engine
def run(parameters,objects,pipeline):
    """preprocess the inputfile with RMA or MAS5
       using preprocess.py will generate a output file"""
    #preprocess the cel file to text signal file
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    import Betsy_config
    PREPROCESS_path = Betsy_config.PREPROCESS
    PREPROCESS_BIN = module_utils.which(PREPROCESS_path)
    assert PREPROCESS_BIN,'cannot find the %s' %PREPROCESS_path
    command = ['python', PREPROCESS_BIN, parameters['preprocess'].upper(), 
               single_object.identifier]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        if not "Loading required package: Biobase" in error_message:
            raise ValueError(error_message)
    outputfiles = os.listdir(os.getcwd())
    for i in outputfiles:
        if i.endswith('.mas5') and not i.endswith('.l2.mas5'):
            outputfile = i
        elif i.endswith('.rma'):
            outputfile = i
    os.rename(outputfile,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for preprocess fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_preprocess_' + original_file + '.jeffs'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'cel_files','contents',['v3_4'])
    assert os.path.exists(single_object.identifier),(
        'the input file %s for preprocess does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
