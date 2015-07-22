#convert_family_soft_to_rename.py
import os
from Betsy import module_utils, userfile
from genomicode import config, arrayplatformlib
from Betsy import bie3
from Betsy import rulebase
from genomicode.filelib import openfh


def extract_sample2desc(GSEID, filename):
    title_dict = {}
    description_dict = {}
    id = None
    for line in openfh(filename):
        if line.startswith("^SAMPLE"):
            assert id is None, "problem with %s" % filename
            id = line.strip().split()[2]
        elif line.startswith("!Sample_description"):
            assert id is not None, "problem with %s" % filename
            title = line.strip().split(None, 2)[2]
            x = id, title
            description_dict[id] = title
        elif line.startswith("!Sample_title"):
            assert id is not None, "problem with %s" % filename
            title = line.strip().split(None, 2)[2]
            x = id, "Title: %s" % title
            title_dict[id] = title
        elif line.startswith("!sample_table_end"):
            id = None
    return title_dict, description_dict


def convert_family_to_relabel_file(relabel_dict, outfile):
    f = file(outfile, 'w')
    f.write('\t'.join(['Name', 'NewName']) + '\n')
    for key in relabel_dict:
        f.write('\t'.join([key, relabel_dict[key]]) + '\n')
    f.close()


def run(network, antecedents, out_attributes, user_options, num_cores):
    """convert family soft file to RenameFile"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    GSEID = user_options['GSEID']
    title, description = extract_sample2desc(GSEID, in_data.identifier)
    if out_attributes['labels_from'] == 'title':
        convert_family_to_relabel_file(title, outfile)
    else:
        convert_family_to_relabel_file(description, outfile)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_family_soft_to_rename does not exists'
        % outfile)
    out_node = bie3.Data(rulebase.RenameFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'rename_file_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)

    return data_node
