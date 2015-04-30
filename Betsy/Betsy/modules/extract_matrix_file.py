#extract_matrix_file.py
import os
import shutil
from genomicode import affyio
from ftplib import FTP
from Betsy import module_utils, bie3, rulebase
import string
import gzip


def run(data_node, parameters, user_input, network, num_cores):
    """extract the matrix file from the expression files"""
    from genomicode import affyio
    outfile = name_outfile(data_node, user_input)
    directory = module_utils.unzip_if_zip(data_node.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for filename in filenames:
        if 'series_matrix.txt' in filename:
            fileloc = os.path.join(directory, filename)
            outname = os.path.join(outfile, filename)
            extract_expression_file(fileloc, outname)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for extract_matrix_file fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'SignalFile_Postprocess_' + original_file
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


def extract_expression_file(fileloc, outname):
    """extract the matrix from the series_matrix file"""
    text_list = open(fileloc, 'r').readlines()
    start = text_list.index('!series_matrix_table_begin')
    end = test_list.index('!series_matrix_table_end')
    f = file(outname, 'w')
    for i in range(start + 1, end):
        line = text_list[i].replace('"', '')
        f.write(line)
    f.close()
