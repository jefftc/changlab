#extract_CEL_files.py
import os
import shutil
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(data_node, parameters, user_input, network, num_cores):
    """extract the cel files with cc or v3_4"""
    from genomicode import affyio
    outfile = name_outfile(data_node, user_input)
    directory = module_utils.unzip_if_zip(data_node.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
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
            'the output file %s for extract_CEL_files fails' % outfile
        )
        out_node = bie3.Data(rulebase.CELFiles, **parameters)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    else:
        assert ValueError('There is no cel file in the input.')


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'cel_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_node):
    return parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node
