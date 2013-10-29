#add_crossplatform_probeid.py

import os
import shutil
import subprocess
import arrayio
from Betsy import module_utils
from genomicode import jmath, Matrix, arrayplatformlib, config
from Betsy import bie,rulebase
 

def run(data_node,parameters,network):
    outfile = name_outfile(data_node)
    DATA = arrayio.read(data_node.attributes['filename'])
    chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
    platform = parameters['platform']
    assert arrayplatformlib.get_bm_attribute(platform), (
        'the desire platform %s is not recognized by Betsy' % platform)
    if chipname == platform:
        shutil.copyfile(data_node.attributes['filename'], outfile)
    else:
        Annot_path = config.annotate_matrix
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN, 'cannot find the %s' % Annot_path
        command = ['python', Annot_BIN, '-f', data_node.attributes['filename'],
                   '-o', outfile, "--platform", platform]
        process = subprocess.Popen(command, shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for add_crossplatform_probeid fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)





