#mark_duplcates.py
import os
from Betsy import module_utils,bie3, rulebase
import subprocess
from genomicode import config

def run(data_node, parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    mark_duplicates_path = config.Mark_duplicates
    assert os.path.exists(mark_duplicates_path),'cannot find the %s' %mark_duplicates_path
    command = ['java','-Xmx5g','-jar',mark_duplicates_path,'I='+data_node.identifier,
                'O='+outfile,
                'METRICS_FILE=metricsFile',
                 'VALIDATION_STRINGENCY=LENIENT', 
                'REMOVE_DUPLICATES=true']
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    process.wait()
    #error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for mark_duplcates does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.BamFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def get_out_attributes(parameters,data_object):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'marked_duplicates_' + original_file+ '.bam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='BamFile')
    
    return data_node

