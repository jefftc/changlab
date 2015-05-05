#rank_genes_by_class_neighbors.py

import os
from genomicode import config
import subprocess
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    import arrayio
    tmp = os.path.join(os.getcwd(), 'tmp.txt')
    f = file(tmp, 'w')
    M = arrayio.read(data_node.identifier)
    M_c = arrayio.convert(M, to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c, f)
    f.close()
    module_name = 'ClassNeighbors'
    gp_parameters = dict()
    gp_parameters['data.filename'] = tmp
    gp_parameters['class.filename'] = cls_node.identifier
    if 'cn_num_neighbors' in user_options:
        gp_parameters['num.neighbors'] = str(user_options['cn_num_neighbors'])
    if 'cn_num_perm' in user_options:
        if user_options['cn_num_perm'].isdigit():
            gp_parameters['num.permutations'] = str(user_options['cn_num_perm'])
    if 'cn_user_pval' in user_options:
        if module_utils.is_number(user_options['cn_user_pval']):
            gp_parameters['user.pval'] = str(user_options['cn_user_pval'])

    mean_median = {'mean': '', 'median': '-d'}
    if out_attributes['cn_mean_or_median'] in ['mean', 'median']:
        gp_parameters['mean.or.median'
                     ] = mean_median[out_attributes['cn_mean_or_median']]

    p = {'t_test': '', 'snr': '-S'}
    if out_attributes['cn_ttest_or_snr'] in p.values():
        gp_parameters['ttest.or.snr'] = p[out_attributes['cn_ttest_or_snr']]

    if out_attributes['cn_filter_data'] in ['yes', 'no']:
        gp_parameters['filter.data'] = str(out_attributes['cn_filter_data'])
    if 'cn_abs_diff' in user_options:
        if module_utils.is_number(user_options['cn_abs_diff']):
            gp_parameters['min.abs.diff'] = str(user_options['cn_abs_diff'])
    if 'cn_min_threshold' in user_options:
        if module_utils.is_number(user_options['cn_min_threshold']):
            gp_parameters['min.threshold'] = str(
                user_options['cn_min_threshold'])
    if 'cn_max_threshold' in user_options:
        if module_utils.is_number(user_options['cn_max_threshold']):
            gp_parameters['max.threshold'] = str(
                user_options['cn_max_threshold'])
    if 'cn_min_folddiff' in user_options:
        if module_utils.is_number(user_options['cn_min_folddiff']):
            gp_parameters['min.fold.diff'] = str(
                user_options['cn_min_folddiff'])

    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module, 'cannot find the %s' % gp_path
    download_directory = os.path.join(os.getcwd(), 'class_neighbors_result')
    command = [gp_module, module_name, '-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters', key + ':' + gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory), (
        'there is no output directory for class_neighbors'
    )
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files, 'gene_pattern get error'
    os.remove(tmp)
    gene_list = []
    for result_file in result_files:
        if result_file.endswith('.odf'):
            f = file(os.path.join(download_directory, result_file), 'r')
            text = f.read()
            text = text.split('\n')
            f.close()
            numline = 8
            startline = 14
            assert text[numline].startswith(
                'NumNeighbors'), 'the odf file format is not right'
            number_gene = int(text[numline].split('=')[1])
            assert text[startline].startswith('1'), 'the start line is not right'

            for line in text[startline:startline + number_gene]:
                lines = line.split('\t')
                gene_list.append(lines[10])
    f = file(outfile, 'w')
    f.write('\t'.join(gene_list))
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for rank_genes_by_class_neighbors fails' % outfile
    )
    out_node = bie3.Data(rulebase.GeneListFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='_SignalFile_Order')
    filter2 = module_utils.AntecedentFilter(datatype_name='ClassLabelFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'gene_list' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
