#convert_signal_to_gct.py
import os
from Betsy import module_utils
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """convert signal file to gct format"""
    import arrayio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    f = file(outfile, 'w')
    M = arrayio.read(in_data.identifier)
    M_c = arrayio.convert(M, to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signal_to_gct fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Filter, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_' + original_file + '.gct'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):

    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
