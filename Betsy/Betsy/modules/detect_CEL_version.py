#detect_CEL_version.py
import os
from Betsy import module_utils
import shutil
from genomicode import affyio
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """convert the cel file with ccl or v3_4 to v3_4"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    new_parameters = get_out_attributes(out_attributes, in_data)
    shutil.copytree(in_data.identifier, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for detect_CEL_version' % outfile
    )
    out_node = bie3.Data(rulebase.CELFiles, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'cel_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def set_out_attributes(antecedents, out_attributes):
    filenames = os.listdir(antecedents.identifier)
    ver_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(antecedents.identifier, filename)
            cel_v = affyio.guess_cel_version(fileloc)
            ver_list.append(cel_v)
    new_parameters = out_attributes.copy()
    if 'cc1' in ver_list:
        new_parameters['version'] = 'cc'
    elif 'v3' in ver_list:
        new_parameters['version'] = 'v3_v4'
    elif 'v4' in ver_list:
        new_parameters['version'] = 'v3_v4'
    else:
        raise ValueError('the cel file can only be cc,v3,v4')
    return new_parameters


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
