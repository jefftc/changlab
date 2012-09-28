#calc_diffexp_with_sam.py
import subprocess
import shutil
import os
import arrayio
from genomicode import jmath
from Betsy import module_utils, read_label_file


def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),(
        'cannot find label_file %s'%label_file.identifier)
    label,label_line,second_line = read_label_file.read(label_file.identifier)
    class_num = len(label)
    assert class_num == 2, 'the number of class in %s is not 2'%label_file.identifier
    delta = 0
    foldchange = 0
    if 'sam_delta' in parameters.keys():
        delta = float(parameters['sam_delta'])
    if 'sam_foldchange' in parameters.keys():
        foldchange = float(parameters['sam_foldchange'])
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    sam_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'sam_script.py')
    cmd = ['python', sam_script, single_object.identifier,
           label_file.identifier, outfile, str(delta), str(foldchange)]
    process = subprocess.Popen(cmd, shell=False,
                     stdout = subprocess.PIPE,
                     stderr = subprocess.STDOUT)
    process.wait()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for calc_diffexp_with_sam fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile)
    return new_objects

def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'sam_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters, objects, 'signal_file', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for calc_diffexp_with_sam does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    parameters = module_utils.renew_parameters(parameters, ['status'])
    new_object = module_utils.DataObject(
        'differential_expressed_genes', [parameters['contents']], outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
   
