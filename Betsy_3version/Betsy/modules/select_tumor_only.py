#select_tumor_only.py
import os
import subprocess
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils
from genomicode import config

def extract_files(gzfile):
    assert gzfile.lower().endswith('tar.gz')
    import tarfile
    tfile = tarfile.open(gzfile, 'r:gz')
    gzname = os.path.split(gzfile)[-1]
    newdir = os.path.join(os.getcwd(),gzname[:-7])
    tfile.extractall(newdir)
    folder = os.listdir(newdir)
    directory = os.path.join(newdir,folder[0])
    assert os.path.exists(directory)
    files = os.listdir(directory)
    for filename in files:
        if filename.lower().endswith('.txt') and filename != 'MANIFEST.txt':
            return os.path.join(directory,filename)
    return None 

def run(data_node,parameters, user_input,network,num_cores):
    """select tumor sample only """
    outfile = name_outfile(data_node,user_input)
    infile = data_node.identifier
    if data_node.identifier.endswith('tar.gz'):       
        infile = extract_files(data_node.identifier)
    slice_matrix_BIN = config.slice_matrix
    slice_matrix = module_utils.which(slice_matrix_BIN)
    assert slice_matrix, 'cannot find the %s' % slice_matrix_BIN
    tempfile='temp.txt'
    
    process = subprocess.Popen([slice_matrix_BIN,
                                '--reorder_col_alphabetical',infile,
                                ],shell=False,
                                stdout=file(tempfile,'w'),
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tempfile),(
        'the temp file %s for select_tumor_only fails'%tempfile)
    tempfile1='temp1.txt'
    process = subprocess.Popen([slice_matrix_BIN,
                                '--tcga_solid_tumor_only',tempfile,
                                ],shell=False,
                                stdout=file(tempfile1,'w'),
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tempfile1),(
        'the temp file %s for select_tumor_only fails'%tempfile1)
    tempfile2='temp2.txt'
    process = subprocess.Popen([slice_matrix_BIN,'--tcga_relabel_patient_barcodes',tempfile1,
                                ],shell=False,
                                stdout=file(tempfile2,'w'),
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tempfile2),(
        'the output file %s for select_tumor_only fails'%tempfile2)
    process = subprocess.Popen([slice_matrix_BIN,'--remove_duplicate_cols','--zerofill',tempfile2,
                                ],shell=False,
                                stdout=file(outfile,'w'),
                                stderr=subprocess.PIPE)  
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    
    assert module_utils.exists_nz(outfile),(
        'the output file %s for select_tumor_only fails'%outfile)
    
    os.remove(tempfile)
    os.remove(tempfile1)
    os.remove(tempfile2)
 
    out_node = bie3.Data(rulebase.TCGAFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    original_file = original_file.replace('.rar','')
    filename = 'TCGAFile_' + original_file+'.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters
    

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
