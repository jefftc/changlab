#preprocess_rma.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    """preprocess the inputfile with RMA 
       using preprocess.py will generate a output file"""
    in_data = antecedents
    #preprocess the cel file to text signal file
    outfile = name_outfile(in_data, user_options)
    PREPROCESS_path = config.preprocess
    PREPROCESS_BIN = module_utils.which(PREPROCESS_path)
    assert PREPROCESS_BIN, 'cannot find the %s' % PREPROCESS_path
    command = ['python', PREPROCESS_BIN, 'RMA', in_data.identifier]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        if not "Loading required package: Biobase" in error_message:
            raise ValueError(error_message)
    
    outputfiles = os.listdir(os.getcwd())
    outputfile = None
    for i in outputfiles:
        if i.endswith('.rma'):
            outputfile = i
    
    os.rename(outputfile, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for preprocess_rma fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_rma_' + original_file + '.jeffs'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
