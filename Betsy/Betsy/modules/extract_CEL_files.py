#extract_CEL_files.py
import os
import shutil
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    """extract the cel files with cc or v3_4"""
    from genomicode import affyio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
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
        out_node = bie3.Data(rulebase.CELFiles, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    else:
        assert ValueError('There is no cel file in the input.')



def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'cel_files_' + original_file
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
