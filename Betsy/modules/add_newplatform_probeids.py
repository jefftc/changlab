#add_newplatform_probeids.py

import os
import module_utils
import subprocess
import Betsy_config
from genomicode import arrayannot,jmath,Matrix
import shutil
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    id,chipname = arrayannot.guess_platform(single_object.identifier)
    platform = parameters['platform']
    assert platform in platforms,'the desire platform  %s is not recognized by Betsy'%platform
    if chipname == platform:
        shutil.copyfile(single_object.identifier,outfile)
    else:
        Annot_path = Betsy_config.ANNOTATE_MATRIX
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN,'cannot find the %s' %Annot_path
        command = ['python', Annot_BIN,'-f',single_object.identifier,'-o',outfile,"--platform",
               platform]
        process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for add_newplatform_probeids fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal' + '_u133plus2_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for add_newplatform_probeids'%single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


platforms=['HG_U133_Plus_2','HG_U133B',
           'HG_U133A','HG_U133A_2',
           'HG_U95A','HumanHT_12',
           'HG_U95Av2','entrez_ID_human',
           'entrez_ID_symbol_human','Hu6800',
           'Mouse430A_2', 'MG_U74Cv2',
           'Mu11KsubB','Mu11KsubA',
           'MG_U74Av2','Mouse430_2','MG_U74Bv2',
           'entrez_ID_mouse','MouseRef_8',
           'entrez_ID_symbol_mouse',
           'RG_U34A','RAE230A']

