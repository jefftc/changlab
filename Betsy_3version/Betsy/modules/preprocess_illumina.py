#preprocess_illumina.py

import shutil
import os
import zipfile
import subprocess
import arrayio
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def zip_directory(dir, zip_file):
    zip = zipfile.ZipFile(zip_file, 'w',
                          compression = zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()

def run(data_node, parameters, user_input, network):
    outfile = name_outfile(data_node,user_input)
    module_name = 'IlluminaExpressionFileCreator'
    gp_parameters = dict()
    if zipfile.is_zipfile(data_node.identifier):
        gp_parameters['idat.zip'] = data_node.identifier
    else:
        zipfile_name = os.path.split(data_node.identifier)[-1]+'.zip'
        zip_directory(data_node.identifier,zipfile_name)
        gp_parameters['idat.zip'] = os.path.join(os.getcwd(),zipfile_name)
    gp_parameters['manifest'] = 'HumanHT-12_V4_0_R2_15002873_B.txt'
    gp_parameters['chip'] = 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'
    manifiest = ['HumanHT-12_V3_0_R2_11283641_A.txt',
                    'HumanHT-12_V4_0_R2_15002873_B.txt',
                    'HumanHT-12_V3_0_R3_11283641_A.txt',
                    'HumanHT-12_V4_0_R1_15002873_B.txt',
                    'HumanMI_V1_R2_XS0000122-MAP.txt',
                    'HumanMI_V2_R0_XS0000124-MAP.txt',
                    'HumanRef-8_V2_0_R4_11223162_A.txt',
                    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
                    'HumanRef-8_V3_0_R2_11282963_A.txt',
                    'HumanRef-8_V3_0_R3_11282963_A.txt',
                    'HumanWG-6_V2_0_R4_11223189_A.txt',
                    'HumanWG-6_V3_0_R2_11282955_A.txt',
                    'HumanWG-6_V3_0_R3_11282955_A.txt',
                    'MouseMI_V1_R2_XS0000127-MAP.txt',
                    'MouseMI_V2_R0_XS0000129-MAP.txt',
                    'MouseRef-8_V1_1_R4_11234312_A.txt',
                    'MouseRef-8_V2_0_R2_11278551_A.txt',
                    'MouseRef-8_V2_0_R3_11278551_A.txt',
                    'MouseWG-6_V1_1_R4_11234304_A.txt',
                    'MouseWG-6_V2_0_R2_11278593_A.txt',
                    'MouseWG-6_V2_0_R3_11278593_A.txt',
                    'RatRef-12_V1_0_R5_11222119_A.txt']
    chip = ['ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
                'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
                'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
                'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
                'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
                'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
                'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
                'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
                'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
                'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
                'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
                'ilmn_RatRef_12_V1_0_R5_11222119_A.chip']
    if 'ill_bg_mode' in parameters.keys():
        assert parameters['ill_bg_mode'] in ['true','false'],'ill_bg_mode should be ill_yes or ill_no'
        gp_parameters['background.subtraction.mode'] = parameters['ill_bg_mode']
        
    if 'ill_chip' in parameters.keys():
        assert  gp_parameters['ill_chip'] in chip,'ill_chip is not correct'
        gp_parameters['chip'] = str(parameters['ill_chip'])
        
    if 'ill_manifest' in parameters.keys():
        assert  gp_parameters['ill_manifest'] in manifiest,'ill_manifest is not correct'
        gp_parameters['manifest'] = str(parameters['ill_manifest'])
        
    if 'ill_coll_mode' in parameters.keys():
        assert  gp_parameters['ill_coll_mode'] in ['none','max','median'],'ill_coll_mode is not correct'
        gp_parameters['collapse.mode'] = str(parameters['ill_coll_mode'])
        
    if 'ill_clm' in user_input.keys():
        gp_parameters['clm'] = str(parameters['ill_clm'])
        
    if 'ill_custom_chip' in user_input.keys():
        gp_parameters['chip'] = str(parameters['ill_custom_chip'])

    if 'ill_custom_manifest' in user_input.keys():
        gp_parameters['custom.manifest'] = str(
            parameters['ill_custom_manifest'])
    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module,'cannot find the %s' %gp_path
    download_directory = os.path.join(os.getcwd(),'illumina_result')
    command = [gp_module, module_name, '-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    goal_file = None
    assert os.path.exists(download_directory),(
        'there is no output directory for illumina')
    #result_files = os.listdir(download_directory)
    #assert 'stderr.txt' not in result_files,('gene_pattern get error '+
    #        'The contents of stderr.txt is:'+
    #        file(os.path.join(download_directory,'stderr.txt')).read())
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for result_file in result_files:
        if result_file == 'System.out':
            continue
        result_path = os.path.join(outfile,result_file)
        M = arrayio.read(os.path.join(download_directory,result_file))
        a = M._col_names['_SAMPLE_NAME']
        b = sorted(a)
        index = []
        for i in b:
            index.append(a.index(i))
        M_new = M.matrix(None,index)
        f = file(result_path,'w')
        arrayio.gct_format.write(M_new,f)
        f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for illumina fails'%outfile)
    out_node = bie3.Data(rulebase.ILLUFolder,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'IlluFolder_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
