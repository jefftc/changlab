#extract_illumina_idat_files.py

import shutil
import os
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, network):
    outfile = name_outfile(data_node)
    directory = module_utils.unzip_if_zip(data_node.attributes['filename'])
    illumina_file = []
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for filename in filenames:
        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
            continue
        if filename.endswith('.idat'):
            illumina_file.append(filename)
    if illumina_file:
        os.mkdir(outfile)
        for filename in illumina_file:
            if filename[:-5].endswith('_Grn'):
                newfilename = filename[:-9] + filename[-5:]
            else:
                newfilename = filename
            old_file = os.path.join(directory, filename)
            new_file = os.path.join(outfile, newfilename)
            shutil.copyfile(old_file, new_file)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_illumina_idat_files fails'
            % outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.IDATFiles,**new_parameters)
        return out_node
    else:
        print 'There is no illumina idat file in the input.'
        return None


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'idat_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
