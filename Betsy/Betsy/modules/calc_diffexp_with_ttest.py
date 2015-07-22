#calc_diffexp_with_ttest.py

import os
from Betsy import module_utils, bie3, rulebase
from genomicode import config
import subprocess


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    diffexp_bin = config.find_diffexp_genes
    assert os.path.exists(diffexp_bin)
    cmd = ['python', diffexp_bin, data_node.identifier, '--cls_file',
           cls_node.identifier, '--algorithm', 'ttest']
    if 'diffexp_foldchange_value' in user_options:
        foldchange = float(user_options['diffexp_foldchange_value'])
        cmd = cmd + ['--fold_change', str(foldchange)]
    handle = open(outfile, 'w')
    try:
        process = subprocess.Popen(cmd,
                                   shell=False,
                                   stdout=handle,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    finally:
        handle.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for calc_diffexp_with_ttest fails' % outfile
    )
    out_node = bie3.Data(rulebase.DiffExprFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='SignalFile')
    filter2 = module_utils.AntecedentFilter(datatype_name='ClassLabelFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 't_test_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
