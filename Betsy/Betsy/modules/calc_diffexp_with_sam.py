#calc_diffexp_with_sam.py
import subprocess
import os
from Betsy import module_utils
from Betsy import bie3, rulebase
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    diffexp_bin = config.find_diffexp_genes
    assert os.path.exists(diffexp_bin)
    delta = 1.0
    if 'sam_delta_value' in user_options:
        delta = float(user_options['sam_delta_value'])
    cmd = ['python', diffexp_bin, data_node.identifier, '--cls_file',
           cls_node.identifier, '--algorithm', 'sam', '--sam_delta',
           str(delta)]
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
        'the output file %s for calc_diffexp_with_sam fails' % outfile
    )
    out_node = bie3.Data(rulebase.DiffExprFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='SignalFile')
    cls_node = module_utils.get_identifier(network, module_id, pool,
                                           user_attributes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'sam_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
