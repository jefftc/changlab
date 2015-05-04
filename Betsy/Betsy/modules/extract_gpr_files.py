#extract_gpr_files.py
import os
import shutil
from Betsy import module_utils
from Betsy import gpr_module
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """extract the files that are gpr format"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
    check = []
    newfiles = []
    files = os.listdir(directory)
    assert files, 'The input folder or zip file is empty.'
    for i in files:
        if i == '.DS_Store':
            pass
        else:
            gpr = gpr_module.check_gpr(os.path.join(directory, i))
            check.append(gpr)
            newfiles.append(i)
    
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    
    for i in range(len(check)):
        if check[i]:
            old_file = os.path.join(directory, newfiles[i])
            new_file = os.path.join(outfile, newfiles[i])
            shutil.copyfile(old_file, new_file)
    
    if True in check:
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_gpr fails' % outfile
        )
        out_node = bie3.Data(rulebase.GPRFiles, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    else:
        assert ValueError('There is no gpr file in the input.')



def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'gpr_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
