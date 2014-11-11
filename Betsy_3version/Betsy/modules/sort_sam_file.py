#sort_sam_file.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config

def run(data_node, parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    sortsam_BIN = config.sortsam
    assert os.path.exists(sortsam_BIN), 'cannot find the %s' % sortsam_BIN
    command = ['java','-Xmx5g','-jar',sortsam_BIN,
                'I='+data_node.identifier,'O='+outfile,'SO=coordinate',
               'VALIDATION_STRINGENCY=LENIENT',
                'CREATE_INDEX=true']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for sort_sam_file does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.BamFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)
def get_out_attributes(parameters,data_object):
    return parameters

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'sorted_' + original_file+ '.bam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SamFile')
    
    return data_node




