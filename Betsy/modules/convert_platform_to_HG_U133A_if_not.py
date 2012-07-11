#convert_platform_to_affyu1332_if_not.py

import os
import module_utils
import subprocess
import Betsy_config
from genomicode import arrayannot
import shutil
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    slice_path = Betsy_config.RENAME
    slice_BIN = module_utils.which(slice_path)
    assert slice_BIN,'cannot find the %s' %slice_path
    chipname = arrayannot.guess_chip(single_object.identifier)
    if chipname in ['HumanHT-12','MouseRef-8']:
        mapfile = Betsy_config.MAPPING
        assert os.path.exists(mapfile),'the mapping file %s does not exists'%mapfile
        import arrayio
        M = arrayio.read(single_object.identifier)
        num_header = len(M._row_names.keys())
        m = str(num_header+1)
        command1 = ['python', slice_BIN, '--select_row_annotation',mapfile+',Match,Best for Both',
                    '--select_row_numeric_annotation',mapfile+',Distance,<=1000',
                    '--add_row_annot',mapfile+',Affymetrix.Probe.Set.ID',
                     '--rename_row_annot', 'Affymetrix.Probe.Set.ID,Probe.Set.ID','--move_row_annot',m+',1',
                    single_object.identifier]
        f=file(outfile,'w')
        process = subprocess.Popen(command1,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
                raise ValueError(error_message)
        f.close()
    elif parameters['preprocess'] in ['mas5','rma']:
        shutil.copyfile(single_object.identifier,outfile)
    else:
        raise ValueError('we cannot convert the platform you input to HG_U133A')
    assert module_utils.exists_nz(outfile),'the output file %s\
                                            for convert_platform fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal' + '_affyu1332_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),'the input file %s\
                          for convert_platform'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
