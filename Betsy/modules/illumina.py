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
    

def run(pipeline_parameters,objects,pipeline):
    identifier,single_object = get_identifier(pipeline_parameters,objects)
    outfile,new_objects = get_outfile(pipeline_parameters,objects,pipeline)
    module_name = 'IlluminaExpressionFileCreator'
    parameters = dict()
    if zipfile.is_zipfile(identifier):
        parameters['idat.zip'] = identifier
    else:
        zipfile_name = os.path.split(identifier)[-1]+'.zip'
        zip_directory(identifier,zipfile_name)
        parameters['idat.zip'] = zipfile_name
    parameters['manifest'] = 'HumanHT-12_V4_0_R2_15002873_B.txt'
    parameters['chip'] = 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'
    
    gp_module = Betsy_config.GENEPATTERN
    command = [gp_module, module_name]
    for key in parameters.keys():
        a = ['--parameters',key+':'+ parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    output_text =  process.stdout.read()
    
    out_lines = output_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    goal_file = None
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    for result_file in result_files:
        assert result_file.endswith('.gct')
        if '-controls' in result_file:
            if 'controls' in pipeline_parameters['preprocess']:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
        else:
            if 'controls' not in pipeline_parameters['preprocess']:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
                
    os.rename(goal_file,outfile)
    module_utils.write_Betsy_parameters_file(
                    parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'geo_dataset','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId','signal_file',pipeline)
    
def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'geo_dataset','Contents,DatasetId')
