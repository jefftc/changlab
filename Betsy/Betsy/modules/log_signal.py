#log_signal.py
import os
from Betsy import module_utils
from genomicode import binreg
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """log the input file"""
    import arrayio
    import math
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    assert not binreg.is_logged_array_data(M), 'the file is logged'
    for i in range(len(M._X)):
        for j in range(len(M._X[i])):
            if M._X[i][j] is not None:
                if float(M._X[i][j]) < 1:
                    M._X[i][j] = 1
                M._X[i][j] = math.log(float(M._X[i][j]), 2)
    
    f = file(outfile, 'w')
    M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for log_signal fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_log_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
