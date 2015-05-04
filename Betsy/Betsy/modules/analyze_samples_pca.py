#analyze_samples_pca.py
import os
from genomicode import pcalib
import arrayio
from Betsy import module_utils
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    X = M._X
    N = 500
    if 'pca_gene_num' in user_options:
        N = int(user_options['pca_gene_num'])
    
    N = min(N, M.nrow())
    index = pcalib.select_genes_var(X, N)
    M_new = M.matrix(index, None)
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for analyze_samples_pca fails' % outfile
    )
    out_node = bie3.Data(rulebase.PcaAnalysis, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'pca_matrix_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
