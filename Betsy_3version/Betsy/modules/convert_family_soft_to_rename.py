#convert_family_soft_to_rename.py
import os
from Betsy import module_utils, userfile
from genomicode import config, arrayplatformlib
from Betsy import bie3
from Betsy import rulebase
from genomicode.filelib import openfh


def extract_sample2desc(**params):
    title_dict = {}
    description_dict = {}
    for name in params:
        exec "%s = params['%s']" % (name, name)
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


def run(data_node, parameters, user_input, network, num_cores):
    """convert family soft file to RenameFile"""
    outfile = name_outfile(data_node, user_input)
    GSEID = user_input['GSEID']
    title, description = extract_sample2desc(GSEID=GSEID,
                                             filename=data_node.identifier)
    if parameters['labels_from'] == 'title':
        convert_family_to_relabel_file(title, outfile)
    else:
        convert_family_to_relabel_file(description, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_family_soft_to_rename does not exists'
        % outfile)
    out_node = bie3.Data(rulebase.RenameFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_object, user_input):
    original_file = module_utils.get_inputid(data_object.identifier)
    filename = 'rename_file_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)


def get_out_attributes(parameters, data_object):
    return parameters


def find_antecedents(network, module_id, pool, parameters, user_attributes):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)

    return data_node
