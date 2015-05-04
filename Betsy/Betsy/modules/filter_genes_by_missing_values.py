#filter_genes_by_missing_values.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    import arrayio
    f_out = file(outfile, 'w')
    M = arrayio.read(in_data.identifier)
    I_good = []
    #get the percentage of gene filter
    percent = float(user_options['filter_value']) / 100
    for i in range(M.dim()[0]):
        missing_count = 0
        for j in range(M.dim()[1]):
            if M._X[i][j] in [None, 'NA']:
                missing_count = missing_count + 1
        if float(missing_count) / M.dim()[1] < percent:
            I_good.append(i)
    
    M_c = M.matrix(I_good, None)
    arrayio.tab_delimited_format.write(M_c, f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for gene_filter fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Impute, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_filter_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
