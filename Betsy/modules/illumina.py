#illumina.py
import module_utils
import shutil
import os
import Betsy_config
import zipfile
import subprocess
import plot_biotin
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

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    module_name = 'IlluminaExpressionFileCreator'
    gp_parameters = dict()
    if zipfile.is_zipfile(identifier):
        gp_parameters['idat.zip'] = identifier
    else:
        zipfile_name = os.path.split(identifier)[-1]+'.zip'
        zip_directory(identifier,zipfile_name)
        gp_parameters['idat.zip'] = zipfile_name
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
        assert parameters['ill_bg_mode'] in ['ill_yes','ill_no'],'ill_bg_mode should be ill_yes or ill_no'
        p = {'ill_yes':'true','ill_no':'false'}
        gp_parameters['background.subtraction.mode'] = p[parameters['ill_bg_mode']]
        
    if 'ill_chip' in parameters.keys():
        assert  gp_parameters['ill_chip'] in chip,'ill_chip is not correct'
        gp_parameters['chip'] = str(parameters['ill_chip'])
        
    if 'ill_manifest' in parameters.keys():
        assert  gp_parameters['ill_manifest'] in manifiest,'ill_manifest is not correct'
        gp_parameters['manifest'] = str(parameters['ill_manifest'])
        
    if 'ill_coll_mode' in parameters.keys():
        assert  gp_parameters['ill_coll_mode'] in ['ill_none','ill_max','ill_median'],'ill_coll_mode is not correct'
        gp_parameters['collapse.mode'] = str(parameters['ill_coll_mode'])[4:]
        
    if 'ill_clm' in parameters.keys():
        gp_parameters['clm'] = str(parameters['ill_clm'])
        
    if 'ill_custom_chip' in parameters.keys():
        gp_parameters['chip'] = str(parameters['ill_custom_chip'])

    if 'ill_custom_manifest' in parameters.keys():
        gp_parameters['custom.manifest'] = str(parameters['ill_custom_manifest'])

    gp_path = Betsy_config.GENEPATTERN
    gp_module = module_utils.which(gp_path)
    assert gp_module,'cannot find the %s' %gp_path
    command = [gp_module, module_name]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    out_text =  process.stdout.read()
    out_lines = out_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    goal_file = None
    
    assert os.path.exists(download_directory),'there is no output directory for illumina'
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    for result_file in result_files:
        assert result_file.endswith('.gct')
        if '-controls' in result_file:
            if 'controls' in parameters['preprocess']:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
        else:
            if 'controls' not in parameters['preprocess']:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
                
    os.rename(goal_file,outfile)
    assert module_utils.exists_nz(outfile),'the output file %s\
                         for illumina fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                    parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'idat_files','contents',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'idat_files','contents')
    assert os.path.exists(identifier),'the input file %s \
             for illumina does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
