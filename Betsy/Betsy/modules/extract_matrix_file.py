#extract_matrix_file.py
import os
#from genomicode import affyio
from ftplib import FTP
from Betsy import module_utils, bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """extract the matrix file from the expression files"""
    #from genomicode import affyio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
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
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'SignalFile_Postprocess_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def extract_expression_file(fileloc, outname):
    """extract the matrix from the series_matrix file"""
    text_list = open(fileloc, 'r').readlines()
    start = text_list.index('!series_matrix_table_begin')
    end = text_list.index('!series_matrix_table_end')
    f = file(outname, 'w')
    for i in range(start + 1, end):
        line = text_list[i].replace('"', '')
        f.write(line)
    f.close()
