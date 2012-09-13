#annotate_with_gsea.py
import module_utils
import shutil
import os
from genomicode import jmath,arrayannot
import config
import subprocess
import read_label_file
import arrayio
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),(
        'cannot find label_file %s for gsea'%label_file.identifier)
    module_name = 'GSEA'
    gp_parameters = dict()
    gp_parameters['expression.dataset'] = single_object.identifier
    gp_parameters['phenotype.labels'] = label_file.identifier
    M = arrayio.read(single_object.identifier)
    x = arrayannot.identify_all_platforms_of_matrix(M)
    chipname = x[0][1]
    gp_parameters['chip.platform']= chipname+'.chip'
    gp_path = config.GENEPATTERN
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
    assert os.path.exists(download_directory),'there is no output directory for GSEA'
    gp_files = os.listdir(download_directory)
    assert 'stderr.txt' not in gp_files,'gene_pattern get error'
    for gp_file in gp_files:
        if gp_file.endswith('.zip'):
            shutil.copy(os.path.join(download_directory,gp_file),outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for gsea fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the train file %s for signature_analysis does not exist'
        %single_object.identifier)
    return single_object

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'gsea_'+original_file+'.zip'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signature_analysis',parameters,objects,single_object)
    return new_objects
