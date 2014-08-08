#extract_rna_files.py
import os
from Betsy import module_utils,bie3,rulebase
import shutil
import gzip


def run(data_node,parameters,user_input,network):
    """extract the rna seq files"""
    outfile = name_outfile(data_node,user_input)
    
    directory = module_utils.unzip_if_zip(data_node.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    format_types = ['fa','fastq','sam','bam']
    for format_type in format_types:
        for filename in filenames:
            if filename == '.DS_Store':
                continue
            fileloc = os.path.join(data_node.identifier, filename)
            if fileloc.endswith(format_type+'.gz'):
                newfname = os.path.splitext(filename)[0]
                new_file = module_utils.gunzip(fileloc)
            elif fileloc.endswith(format_type):
                new_file = fileloc
                newfname = filename
                shutil.copyfile(new_file, os.path.join(outfile, newfname))
            if fileloc.endswith('.gz'):
                os.remove(new_file)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for extract_rna_files fails' % outfile)
    out_node = bie3.Data(rulebase.RNA_SeqFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'rna_seq_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    return parameters

def find_antecedents(network, module_id,data_nodes, parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='RNA_SeqFile')
    
    return data_node
