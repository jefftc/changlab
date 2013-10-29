#convert_CEL_to_v3_v4.py
import os
from Betsy import module_utils
import shutil
import gzip
from Betsy import bie, rulebase

def run(data_node,parameters, network):
    """convert the cel file with ccl or v3_4 to v3_4"""
    from genomicode import affyio
    outfile = name_outfile(data_node)
    filenames = os.listdir(data_node.attributes['filename'])
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    ver_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(data_node.attributes['filename'], filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if fileloc.endswith('.gz'):
                newcelfname = os.path.splitext(filename)[0]
                cel_file = module_utils.gunzip(fileloc)
            else:
                cel_file = fileloc
                newcelfname = filename
            if cel_v == 'cc1':
                f = file(os.path.join(outfile, newcelfname), 'w')
                affyio.convert_cel_cc1_to_3(cel_file, f)
                f.close()
            elif cel_v in ('v3', 'v4'):
                shutil.copyfile(cel_file, os.path.join(outfile, newcelfname))
            if fileloc.endswith('.gz'):
                os.remove(cel_file)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_CEL_to_v3_v4 fails' % outfile)
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
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node






