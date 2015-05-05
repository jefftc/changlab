#reorder_genes_diff.py

import os
import arrayio
from genomicode import arrayplatformlib, config
from Betsy import bie3
from Betsy import rulebase
from Betsy import gene_ranking
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, gene_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    #read the gene order list
    gene_list = open(gene_node.identifier, 'r').read().split()
    M = arrayio.read(data_node.identifier)
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    if x:
        id = x[0][0]
        platform = x[0][1]
        chip = arrayplatformlib.identify_platform_of_annotations(gene_list)
        if not chip:
            chip = []
        signal_file = data_node.identifier
        if platform == chip:
            tmpfile = data_node.identifier
        else:
            platform_name = 'unknown_platform'
            if 'platform_name' in user_options:
                platform_name = user_options['platform_name']
            if platform_name in chip:  #, 'unknown_platform':
                import subprocess
                Annot_path = config.annotate_matrix
                Annot_BIN = module_utils.which(Annot_path)
                assert Annot_BIN, 'cannot find the %s' % Annot_path
                signal_file = 'tmp'
                command = [
                    'python',
                    Annot_BIN,
                    # Needs to be tested.
                    #'-f', single_object.identifier,
                    '-f', data_node.identifier,
                    '-o', signal_file,
                    "--platform", chip,
                    ]
                process = subprocess.Popen(command,
                                           shell=False,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                error_message = process.communicate()[1]
                if error_message:
                    raise ValueError(error_message)
                assert module_utils.exists_nz(
                    signal_file), 'the platform conversion fails'
                id = out_attributes['platform']
                M = arrayio.read(signal_file)
            elif platform_name == platform:
                # Needs to be tested.
                infile = gene_node.identifier
                #infile = gene_list_file.identifier
                f = file(infile, 'rU')
                genes = f.readlines()
                f.close()
                gene_list = module_utils.convert_gene_list_platform(genes,
                                                                    platform)
    else:
        id = M._row_order[0]
    original_list = M._row_names[id]
    #get the order index and write to the outout file
    indexlist = gene_ranking.find_sorted_index(original_list, gene_list)
    M_new = M.matrix(indexlist, None)
    f = open(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for reorder_genes fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Order, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='_SignalFile_Order')
    filter2 = module_utils.AntecedentFilter(datatype_name='GeneListFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'gene_reorder' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
