#filter_and_threshold_genes.py
import os
from Betsy import module_utils
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """run preprocessdataset """
    import arrayio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    threshold = 20
    ceiling = 16000
    min_fold_change = 5
    min_delta = 100.0
    M = arrayio.read(in_data.identifier)
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        for j in range(len(X[i])):
            if X[i][j] < threshold:
                M._X[i][j] = threshold
            if X[i][j] > ceiling:
                M._X[i][j] = ceiling
        gene = M._X[i]
        fold_change = max(gene) / float(min(gene))
        delta = max(gene) - min(gene)
        if fold_change >= min_fold_change and delta >= min_delta:
            I_good.append(i)
    
    f = file(outfile, 'w')
    M_c = M.matrix(I_good, None)
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for filter_and_threshold_genes fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_preprocessdataset_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
