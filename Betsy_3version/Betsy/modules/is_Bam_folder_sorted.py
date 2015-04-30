#is_Bam_folder_sorted.py
import os
from Betsy import module_utils
import shutil
import gzip
from Betsy import bie3, rulebase
from genomicode import config


def run(data_node, parameters, user_input, network, num_cores):
    """check bamfiles in bam folder is all sorted"""
    from genomicode import affyio
    outfile = name_outfile(data_node, user_input)
    filenames = os.listdir(data_node.identifier)
    parameters = get_out_attributes(parameters, data_node)
    shutil.copytree(data_node.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for is_Bam_folder_sorted fails' % outfile
    )

    out_node = bie3.Data(rulebase.BamFolder, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'BamFolder_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def get_out_attributes(parameters, data_node):
    new_parameters = parameters.copy()
    filenames = os.listdir(data_node.identifier)
    assert filenames, 'The input folder is empty.'
    sorted_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(data_node.identifier, filename)
            sort_info = get_sort_information(fileloc)
            sorted_list.append(sort_info)
    if False in sorted_list:
        new_parameters['sorted'] = 'no'
    else:
        new_parameters['sorted'] = 'yes'
    return new_parameters


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node


def get_sort_information(filename):
    import subprocess
    samtools_BIN = config.samtools
    assert os.path.exists(samtools_BIN), 'cannot find the %s' % samtools_BIN
    command = [samtools_BIN, 'view', '-H', filename]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    out_message = process.communicate()
    if 'error' in out_message[1]:
        raise ValueError(error_message)
    message = out_message[0]
    header = None
    for line in message.split('\n'):
        if line.startswith('@HD'):
            header = line
            break
    sort = header.split('SO:')[-1]
    if sort == 'coordinate':
        return True
    elif sort == 'unsorted':
        return False
    else:
        raise ValueError('cannot determine sorted or not')
