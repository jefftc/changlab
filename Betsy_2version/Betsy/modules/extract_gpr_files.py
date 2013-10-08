#extract_gpr_files.py
import os
import gzip
import shutil
#from Betsy
import module_utils
from Betsy import gpr_module
import bie
import rulebase

def run(data_node,parameters):
    """extract the files that are gpr format"""
    outfile = name_outfile(data_node)
    directory = module_utils.unzip_if_zip(data_node.attributes['filename'])
    check = []
    newfiles = []
    files = os.listdir(directory)
    assert files, 'The input folder or zip file is empty.'
    for i in files:
        if i == '.DS_Store':
            pass
        else:
            gpr = gpr_module.check_gpr(
                os.path.join(directory, i))
            check.append(gpr)
            newfiles.append(i)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for i in range(len(check)):
        if check[i]:
            old_file = os.path.join(directory, newfiles[i])
            new_file = os.path.join(outfile, newfiles[i])
            shutil.copyfile(old_file, new_file)
    if True in check:
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_gpr fails' % outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.GPRFiles,**new_parameters)
        return out_node

    else:
        print 'There is no gpr file in the input.'
        return None



def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'gpr_files_' + original_file
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
