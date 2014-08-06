#call_variants_mpileup.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config

def run(data_node, parameters, user_input, network):
    outfile = name_outfile(data_node,user_input)
    species = parameters['ref']
    if species == 'hg18':
       ref_file = config.hg18_ref
    elif species == 'hg19':
       ref_file = config.hg19_ref
    elif species == 'dm3':
       ref_file = config.dm3_ref
    elif species == 'mm9':
       ref_file= config.mm9_ref
    assert os.path.exists(ref_file),'the ref file %s does not exsits' %ref_file
    #command = ['samtools','mpileup','-uf',ref,single_object.identifier,'|',
    #           'bcftools','view', '-bvcg','-','>',outfile]
    samtools_BIN = config.samtools
    samtools_module = module_utils.which(samtools_BIN)
    assert os.path.exists(samtools_module), 'cannot find the %s' % samtools_BIN
    command = [samtools_BIN,'mpileup','-uf',ref_file,data_node.identifier]
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
        'the output file %s for call_variants_mpileup does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.VcfFile,**parameters)
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
    filename = 'mpileup_' + original_file+ '.bcf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='BamFile')
    
    return data_node









