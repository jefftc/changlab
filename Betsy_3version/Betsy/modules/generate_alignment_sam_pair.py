#generate_alignment_sam_pair.py
import os
from Betsy import module_utils,bie3, rulebase
import subprocess
from genomicode import config


def run(in_nodes, parameters, user_input, network,num_cores):
    fastq1_node, sai1_node, fastq2_node, sai2_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    species = parameters['ref']
    if species == 'hg18':
       ref_file = config.hg18_ref
    elif species == 'hg19':
       ref_file = config.hg19_ref
    elif species == 'dm3':
       ref_file = config.dm3_ref
    elif species == 'mm9':
       ref_file= config.mm9_ref
    else:
        raise ValueError('cannot handle %s' %species)
    assert os.path.exists(ref_file),'the ref_file %s does not exist' % ref_file
    bwa_BIN = config.bwa
    bwa_module = module_utils.which(bwa_BIN)
    assert bwa_module, 'cannot find the %s' % bwa_BIN
    command = [bwa_BIN,'sampe',ref_file, sai1_node.identifier,  sai2_node.identifier,
               fastq1_node.identifier,fastq1_node.identifier]
    f=file(outfile,'w')
    try:
        process=subprocess.Popen(command,shell=False,
                             stdout=f,
                             stderr=subprocess.PIPE)
    finally:
        f.close()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for generate_alignment_sam_pair does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.SamFile,**parameters)
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
    filename = 'generate_alignment_sam' + original_file+ '.sam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    fastq1_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='FastqFile',
                                              optional_key='read',optional_value='pair1')
    sai1_node = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='SaiFile',
                                            optional_key='read',optional_value='pair1')
    fastq2_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='FastqFile',
                                              optional_key='read',optional_value='pair2')
    sai2_node = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='SaiFile',
                                            optional_key='read',optional_value='pair2')
    return fastq1_node, sai1_node, fastq2_node, sai2_node



