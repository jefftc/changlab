#extract_gpr_files.py
import os
import gzip
import shutil
from Betsy import module_utils
from Betsy import gpr_module
from Betsy import bie3,rulebase

def run(data_node,parameters, user_input, network):
    """extract the files that are gpr format"""
    outfile = name_outfile(data_node,user_input)
    directory = module_utils.unzip_if_zip(data_node.identifier)
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
	out_node = bie3.Data(rulebase.GPRFiles,**parameters)
    	out_object = module_utils.DataObject(out_node,outfile)
    	return out_object
    else:
        print 'There is no gpr file in the input.'
        return None



def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'gpr_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
