#normalize_samples_with_quantile.py

import os
from genomicode import quantnorm
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    Y = quantnorm.normalize(M)
    f = file(outfile, 'w')
    Y_c = arrayio.convert(Y, to_format=arrayio.pcl_format)
    arrayio.tab_delimited_format.write(Y_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for quantile fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Merge, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_quantile_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def set_out_attributes(antecedents, out_attributes):
    new_parameters = antecedents.data.attributes.copy()
    new_parameters['quantile_norm'] = 'yes'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
