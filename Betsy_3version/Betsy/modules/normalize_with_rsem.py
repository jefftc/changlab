#normalize_with_rsem.py
import os
from Betsy import module_utils,bie3, rulebase
import subprocess
from genomicode import config
import tempfile
import shutil

def guess_version(chromosomes):
    ref_path = config.chrom2genome
    filenames = os.listdir(ref_path)
    for filename in filenames:
        flag=True
        f=file(os.path.join(ref_path,filename))
        version = f.readlines()
        version = [i.strip() for i in version]
        f.close()
        for chromosome in chromosomes:
            if chromosome not in version:
                flag=False
                break    
        if flag:
            version_name=os.path.splitext(filename)[0]
            version_file=os.path.join(ref_path,os.path.splitext(filename)[0])
            return version_name
    return None
        
def guess_format_and_version(input_file):
    command='samtools view '+ input_file +'|cut -f 3 |uniq'
    text=subprocess.check_output(command,stderr=subprocess.PIPE,shell=True)
    text= text.split()
    text=[i for i in text if i!='*']
    version_name = guess_version(text)
    format_type = os.path.splitext(input_file)[-1]
    return format_type, version_name

        
def preprocess_single_sample(input_file,temp_file):
    #still need to handle format and ref
    format_type,ref = guess_format_and_version(input_file)
    options = ''
    if format_type =='.sam':
        options = '--sam'
    elif format_type == '.bam':
        options = '--bam'
    if ref in ['Ensembl.human','Broad.hg18']:
        ref_file = config.rna_hum
    elif ref == 'Ensembl.mouse':
        ref_file = config.rna_mouse
   # elif ref == 'dm3':
   #     ref_file = config.rna_fly
    else:
        raise ValueError("we cannot handle %s" % ref)
    filename = os.path.split(input_file)[-1]
    ID = os.path.splitext(filename)[0]
    if options:
        command = ['rsem-calculate-expression',options,
                    input_file,ref_file,'--no-bam-output','-p','8',ID]
    else:
##        command = ['rsem-calculate-expression','--paired-end',input_file,
##                   '/home/xchen/NGS/try_RSEM/Priyatansh_048_data/RNA-10_L1_2.fastq',ref_file,
##                    '--no-bam-output','-p','8',ID]
        command = ['rsem-calculate-expression',input_file,
                   ref_file,
                    '--no-bam-output','-p','8',ID]
    process=subprocess.Popen(command,shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    process.wait()
    slice_BIN = config.slice_matrix
    command = ['python',slice_BIN,ID+'.genes.results',
               '--select_col_ids','transcript_id,gene_id,TPM',
               '--replace_col_ids','TPM,'+ ID]
    f=file(temp_file,'w')
    try:
        process=subprocess.Popen(command,shell=False,
                             stdout=f,
                             stderr=subprocess.PIPE)
        process.wait()
    finally:
        f.close()
    error_message = process.communicate()[1]
    if 'error' in error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(temp_file), (
        'the output file %s does not exist'
        % temp_file)

def preprocess_multiple_sample(folder,outfile):
    filenames = os.listdir(folder)
    file_list = []
    for filename in filenames:
        temp_file = tempfile.mkstemp()[1]
        preprocess_single_sample(os.path.join(folder,filename),
                                 temp_file)
        file_list.append(temp_file)
    result_file = file_list[0]
    tmp_list=file_list[:]
    try:
        for filename in file_list[1:]:
            tmp_result=tempfile.mkstemp()[1]
            f=file(tmp_result,'a+')
            try:
                module_utils.merge_two_files(result_file,filename,f)
            finally:
                f.close()
                tmp_list.append(tmp_result)
            result_file = tmp_result
        os.rename(result_file,outfile)
    finally:
        for filename in tmp_list:
            if os.path.exists(filename):
                os.remove(filename)
    
def run(data_node,parameters,user_input,network):
    outfile = name_outfile(data_node,user_input)
    preprocess_multiple_sample(data_node.identifier,outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for normalize_with_rsem does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.SignalFile_Postprocess,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = original_file +'.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes, parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='RNA_SeqFile')
    
    return data_node


