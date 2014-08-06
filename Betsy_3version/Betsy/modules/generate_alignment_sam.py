#generate_alignment_sam.py
import os
from Betsy import module_utils,bie3, rulebase
import subprocess
from genomicode import config


def run(in_nodes, parameters, user_input, network):
    fastq_node,sai_node = in_nodes
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
        raise ValueError('cannot handle %s'%species)
    assert os.path.exists(ref_file),'the ref_file %s does not exist' % ref_file
    bwa_BIN = config.bwa
    bwa_module = module_utils.which(bwa_BIN)
    assert bwa_module, 'cannot find the %s' % bwa_BIN
    command = [bwa_BIN,'samse',ref_file, sai_node.identifier, fastq_node.identifier]
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
        'the output file %s for generate_alignment_sam does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.SamFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    


def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    fastq_node, sai_node = in_nodes
    identifier = sai_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)
def get_out_attributes(parameters,data_object):
    return parameters

def name_outfile(in_nodes,user_input):
    fastq_node, sai_node = in_nodes
    original_file = module_utils.get_inputid(sai_node.identifier)
    filename = 'generate_alignment_sam' + original_file+ '.sam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id,data_nodes,parameters):
    fastq_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='FastqFile')
    sai_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='SaiFile')
    return fastq_node, sai_node



