#run_RNA_SeQC.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config

def run(data_node, parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    directory = module_utils.unzip_if_zip(data_node.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    RNA_SeQC_path = config.RNASeQC
    REF = user_input['RNA_ref']
    GTF = user_input['RNA_gtf']
    assert os.path.exists(RNA_SeQC_path),('cannot find the %s'
                                        %RNA_SeQC_path)
    assert os.path.exists(REF),('cannot find the %s'
                                        %REF)
    assert os.path.exists(GTF),('cannot find the %s'
                                        %GTF)
    for filename in filenames:
        infile = os.path.join(directory, filename)
        outname = os.path.splitext(filename)[0]
        
        command = ['java','-jar',RNA_SeQC_path,'-o',outfile,
               '-r',REF, '-s',outname+'|'+infile+'|NA','-t',GTF]
        process=subprocess.Popen(command,shell=False,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()
        # guess if there is error message in the output
        if 'error' in error_message[1]:
            raise ValueError(error_message)
        outname = os.path.join(outfile,outname)
        assert module_utils.exists_nz(outname), (
            'the output file %s for run_RNA_SeQC does not exist'
            % outname)
    out_node = bie3.Data(rulebase.RNASeQCFile,**parameters)
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
    filename = 'bamFolder_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='BamFolder')
    
    return data_node




