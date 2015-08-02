#is_sam_folder.py
import os
from Betsy import module_utils, bie3, rulebase
import shutil


def run(network, antecedents, out_attributes, user_options, num_cores):
    """check is sam folder"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    
    format_type = 'sam'
    for filename in filenames:
        if filename == '.DS_Store':
            continue
        fileloc = os.path.join(in_data.identifier, filename)
        if fileloc.endswith(format_type + '.gz'):
            newfname = os.path.splitext(filename)[0]
            new_file = module_utils.gunzip(fileloc)
        elif fileloc.endswith(format_type):
            new_file = fileloc
            newfname = filename
            shutil.copyfile(new_file, os.path.join(outfile, newfname))
        if fileloc.endswith('.gz'):
            os.remove(new_file)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for is_sam_folder fails' % outfile
    )
    out_node = bie3.Data(rulebase.RNASeqFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'sam_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    in_data = antecedents
    directory = module_utils.unzip_if_zip(in_data.identifier)
    filenames = os.listdir(directory)
    if directory != antecedents.identifier:
        shutil.rmtree(directory)
    assert filenames, 'The input folder or zip file is empty.'
    format_type = 'sam'
    flag = []
    for filename in filenames:
        if filename == '.DS_Store':
            continue
        fileloc = os.path.join(in_data.identifier, filename)
        if fileloc.endswith(format_type + '.gz'):
            flag.append(True)
        elif fileloc.endswith(format_type):
            flag.append(True)
        else:
            flag.append(False)
    if True in flag:
        out_attributes['format_type'] = 'samfolder'
        return out_attributes
    out_attributes['format_type'] = 'not_samfolder'
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='RNASeqFile')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node