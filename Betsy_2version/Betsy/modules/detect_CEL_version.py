#detect_CEL_version.py
import os
from Betsy import module_utils
import shutil
import gzip
from genomicode import affyio
from Betsy import bie, rulebase

def run(data_node, parameters, network):
    """convert the cel file with ccl or v3_4 to v3_4"""
    outfile = name_outfile(data_node)
    new_parameters = get_out_attributes(parameters,data_node)
    shutil.copytree(data_node.attributes['filename'],outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for detect_CEL_version' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.CELFiles,**new_parameters)
    return out_node
        


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'cel_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def get_out_attributes(parameters,data_node):
    filenames = os.listdir(data_node.attributes['filename'])
    ver_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(data_node.attributes['filename'], filename)
            cel_v = affyio.guess_cel_version(fileloc)
            ver_list.append(cel_v)
    new_parameters = parameters.copy()
    if 'cc1' in ver_list:
        new_parameters['version'] = 'cc'
    elif 'v3' in ver_list:
        new_parameters['version'] = 'v3'
    elif 'v4' in ver_list:
        new_parameters['version'] = 'v4'
    else:
        raise ValueError('the cel file can only be cc,v3,v4')
    return new_parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

