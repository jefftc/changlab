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
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    import Betsy_config
    PREPROCESS_BIN = Betsy_config.PREPROCESS
    command = ['python', PREPROCESS_BIN, parameters['preprocess'].upper(), 
               identifier]
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
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'geo_dataset','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
