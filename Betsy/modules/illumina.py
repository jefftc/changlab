#illumina.py
import module_utils
import shutil
import os
import Betsy_config
import zipfile

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
    

def run(pipeline_parameters,objects):
    identifier,single_object = get_identifier(pipeline_parameters,objects)
    outfile,new_objects = get_outfile(pipeline_parameters,objects)
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
    download_directory = module_utils.run_gp_module(
        module_name,parameters)
    goal_file = None
    result_files = os.listdir(download_directory)
    for result_file in result_files:
        assert result_file.endswith('.gct')
        if 'controls' in pipeline_parameters['preprocess']:
            if '-controls' in result_file:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
        else:
            if '-controls' not in result_file:
                goal_file = os.path.realpath(
                        download_directory+'/'+result_file)
        os.rename(goal_file,outfile)
        module_utils.write_Betsy_parameters_file(
                    parameters,single_object)
        return new_objects
    else:
        return None

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(
        parameters,objects,'geo_dataset','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'geo_dataset','Contents,DatasetId')
