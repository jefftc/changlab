#extract_CEL_files.py
import os
#from Betsy
import module_utils
import shutil
from time import strftime,localtime
import bie
import rulebase

def run(data_node, parameters):
    """extract the cel files with cc or v3_4"""
    from genomicode import affyio
    outfile = name_outfile(data_node)
    directory = module_utils.unzip_if_zip(data_node.attributes['filename'])
    filenames = os.listdir(directory)
    assert filenames,'The input folder or zip file is empty.'
    ver_list = []
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(directory, filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if cel_v in ['cc1', 'v3', 'v4']:
                shutil.copyfile(fileloc, os.path.join(outfile, filename))
                ver_list.append(True)
            else:
                ver_list.append(False)
   
    if True in ver_list:
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_CEL_files fails' % outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.CELFiles,**new_parameters)
        return out_node
    else:
        print 'There is no cel file in the input.'
        return None



def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'cel_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
