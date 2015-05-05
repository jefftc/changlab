#fill_missing_with_zeros.py
import os
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    assert module_utils.is_missing(in_data.identifier), 'no missing values'
    M = arrayio.read(in_data.identifier)
    f_out = file(outfile, 'w')
    for i in range(M.dim()[0]):
        for j in range(M.dim()[1]):
            if M._X[i][j] is None:
                M._X[i][j] = '0'
    
    arrayio.tab_delimited_format.write(M, f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for zero_fill_if_missing fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Impute, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_zero_fill_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
