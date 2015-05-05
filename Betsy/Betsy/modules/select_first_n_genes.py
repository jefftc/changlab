# select_first_n_genes.py
import os
from Betsy import module_utils
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """select a num of genes according to num_features"""
    import arrayio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    f = file(outfile, 'w')
    num_features = 500
    if 'num_features_value' in user_options:
        num_features = int(user_options['num_features_value'])
    
    assert num_features > 0, 'the num_features should be >0'
    M = arrayio.read(in_data.identifier)
    M_c = M.matrix(range(0, num_features), None)
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for select_first_n_genes fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Filter, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_select_n_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):

    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
