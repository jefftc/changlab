#add_U133A_probeid.py

import os
import shutil
import subprocess
import arrayio
from Betsy import module_utils
from genomicode import jmath, Matrix, arrayplatformlib, config
from Betsy import bie3, rulebase


def run(data_node, parameters, user_input, network, num_cores):
    outfile = name_outfile(data_node, user_input)
    DATA = arrayio.read(data_node.identifier)
    chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
    platform = 'HG_U133A'
    assert arrayplatformlib.get_bm_attribute(platform), (
        'the desire platform %s is not recognized by Betsy' % platform
    )
    if chipname == platform:
        shutil.copyfile(data_node.identifier, outfile)
    else:
        Annot_path = config.annotate_matrix
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN, 'cannot find the %s' % Annot_path
        command = ['python', Annot_BIN, data_node.identifier, "--platform",
                   platform, '--min_match_score', '0.89']
        f = file(outfile, 'w')
        try:
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=f,
                                       stderr=subprocess.PIPE)
        finally:
            f.close()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    #change the HG_U133A to the first column
    f = file(outfile, 'r')
    txt = f.readlines()
    f.close()
    header = txt[0].split('\t')
    index = header.index('HG_U133A')
    f = file(outfile, 'w')
    for line in txt:
        line = line.split('\t')
        newline = [line[index]] + line[0:index] + line[index + 1:]
        f.write('\t'.join(newline))
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for add_U133A_probeid fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Annotate, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_node):
    return parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)
