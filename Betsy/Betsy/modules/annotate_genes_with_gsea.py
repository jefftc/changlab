#annotate_genes_with_gsea.py
import subprocess
import shutil
import os
import arrayio
from genomicode import arrayplatformlib, config
from Betsy import module_utils, read_label_file
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    gsea_path = config.gsea
    gsea_module = module_utils.which(gsea_path)
    assert gsea_module, 'cannot find the %s' % gsea_path
    M = arrayio.read(data_node.identifier)
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    chipname = x[0][1]
    assert chipname in arrayplatformlib.platform_to_GSEA_chipname, (
        'we cannot find chipname %s in gsea' % chipname
    )
    chipname_in_gsea = arrayplatformlib.platform_to_GSEA_chipname[chipname]
    platform = chipname_in_gsea + '.chip'
    download_directory = os.path.join(os.getcwd(), 'gsea_result')
    command = [gsea_module, data_node.identifier, '--cls_file',
               cls_node.identifier, '--platform', platform, '--database',
               out_attributes['geneset_database'], '--permutation_type',
               out_attributes['permutation_type'], download_directory]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(download_directory), (
        'there is no output directory for GSEA'
    )
    shutil.copytree(download_directory, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annotate_genes_with_gsea fails' % outfile
    )
    out_node = bie3.Data(rulebase.GseaFile, **out_attributes)
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
    filename = 'gsea_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
