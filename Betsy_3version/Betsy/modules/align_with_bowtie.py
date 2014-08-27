#align_with_bowtie.py
import os
from Betsy import module_utils,bie3, rulebase
import subprocess
from genomicode import config
import tempfile
import shutil



def preprocess_multiple_sample(folder, group_dict, outfile,ref):
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    if ref == 'human':
        ref_file = config.rna_hum
    elif ref == 'mouse':
        ref_file = config.rna_mouse
    else:
        raise ValueError("we cannot handle %s" % ref)
    for sample in group_dict:
        files = group_dict[sample]
        if len(files)==1:
            input_file = os.path.join(folder,files[0])
            command = ['bowtie', '-q','--phred33-quals','-n','2','-e','99999999',
                       '-l','25', '-p', '8' ,'-a', '-m',
                       '200' ,'-S',ref_file,input_file]
        elif len(files)==2:
            input_file1 = os.path.join(folder,files[0])
            input_file2 = os.path.join(folder,files[1])
            command = ['bowtie', '-q','--phred33-quals','-n','2','-e','99999999',
                       '-l','25','-I', '1', '-X' ,'1000', '-p', '8' ,'-a', '-m',
                       '200' ,'-S',ref_file,'-1',input_file1,'-2',input_file2]
        else:
            raise ValueError('number files is not correct')
        outfilename = os.path.join(outfile,sample+'.sam')
        f=file(outfilename,'w')
        try:
            process=subprocess.Popen(command,shell=False,
                                 stdout=f,
                                 stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()[1]
            if 'error' in error_message:
                raise ValueError(error_message)       
        finally:
            f.close()


def process_group_info(group_file):
    f = file(group_file,'r')
    text = f.readlines()
    f.close()
    group_dict = {}
    text = [line.strip() for line in text if line.strip()]
    for line in text:
        words = line.split('\t')
        if len(words)==2: 
            group_dict[words[0]] = [words[1]]
        elif len(words)==3:
            group_dict[words[0]] = [words[2],words[3]]
        else:
            raise ValueError('group file is invalid')
    return group_dict


def run(in_nodes,parameters,user_input,network):
    data_node, group_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    group_dict = process_group_info(group_node.identifier)
    preprocess_multiple_sample(data_node.identifier, group_dict,
                               outfile, parameters['ref'])
    assert module_utils.exists_nz(outfile), (
        'the output file %s for align_with_bowtie does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.SamFolder,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

        
def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,group_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters,user_input)


def name_outfile(in_nodes,user_input):
    data_node,group_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'Samfolder_'+original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,in_nodes):
    return parameters

def find_antecedents(network, module_id,data_nodes, parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='FastqFolder')
    group_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SampleGroupFile')
    
    return data_node, group_node

