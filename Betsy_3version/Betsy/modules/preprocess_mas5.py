#preprocess_mas5.py
import os
from Betsy import module_utils, bie3, rulebase
from genomicode import config
import subprocess


def run(data_node, parameters, user_input, network,num_cores):
    """preprocess the inputfile with  MAS5
       using preprocess.py will generate a output file"""
    #preprocess the cel file to text signal file
    outfile = name_outfile(data_node,user_input)
    PREPROCESS_path = config.preprocess
    PREPROCESS_BIN = module_utils.which(PREPROCESS_path)
    assert PREPROCESS_BIN,'cannot find the %s' %PREPROCESS_path
    command = ['python', PREPROCESS_BIN, 'MAS5', 
               data_node.identifier]
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
    os.rename(outputfile,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for preprocess_mas5 fails'%outfile)
    out_node = bie3.Data(rulebase._SignalFile_Postprocess,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

    
def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_mas5_' + original_file + '.jeffs'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
