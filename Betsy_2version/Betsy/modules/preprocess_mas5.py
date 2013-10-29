#preprocess_mas5.py
import os
from Betsy import module_utils, bie, rulebase
from genomicode import config
import subprocess


def run(data_node, parameters, network):
    """preprocess the inputfile with  MAS5
       using preprocess.py will generate a output file"""
    #preprocess the cel file to text signal file
    outfile = name_outfile(data_node)
    PREPROCESS_path = config.preprocess
    PREPROCESS_BIN = module_utils.which(PREPROCESS_path)
    assert PREPROCESS_BIN,'cannot find the %s' %PREPROCESS_path
    command = ['python', PREPROCESS_BIN, 'MAS5', 
               data_node.attributes['filename']]
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
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node

    
def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_rma_' + original_file + '.jeffs'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
