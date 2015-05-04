#convert_CEL_to_v3_v4.py
import os
from Betsy import module_utils
import shutil
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """convert the cel file with ccl or v3_4 to v3_4"""
    from genomicode import affyio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    filenames = os.listdir(in_data.identifier)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    
    #ver_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(in_data.identifier, filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if fileloc.endswith('.gz'):
                newcelfname = os.path.splitext(filename)[0]
                cel_file = module_utils.gunzip(fileloc)
            else:
                cel_file = fileloc
                newcelfname = filename
            if cel_v == 'cc1':
                f = file(os.path.join(outfile, newcelfname), 'w')
                affyio.convert_cel_cc1_to_3(cel_file, f)
                f.close()
            elif cel_v in ('v3', 'v4'):
                shutil.copyfile(cel_file, os.path.join(outfile, newcelfname))
            if fileloc.endswith('.gz'):
                os.remove(cel_file)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_CEL_to_v3 fails' % outfile
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


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
