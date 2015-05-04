#convert_signalfile_preprocess.py
import os
import shutil
from Betsy import module_utils
from genomicode import binreg
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    from genomicode import filelib
    
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)

    # Copy the file objects so that gzip'd files get properly uncompressed.
    fsrc = filelib.openfh(in_data.identifier)
    fdst = open(outfile, 'w')
    shutil.copyfileobj(fsrc, fdst)
    fdst.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signalfile_preprocess fails' % outfile
    )
    out_node = bie3.Data(rulebase.SignalFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    format_type = antecedents.data.attributes['format']
    filename = 'signal_file_' + original_file + '.' + format_type
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
