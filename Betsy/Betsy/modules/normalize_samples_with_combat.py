#normalize_samples_with_combat.py
import os
import subprocess
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils, read_label_file


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, cls_node = antecedents
    if data_node and cls_node:
        outfile = name_outfile(antecedents, user_options)
        result, label_line, second_line = read_label_file.read(
            cls_node.identifier)
        assert len(
            result) >= 2, 'for combat,there should be equal or larger than 2 classes'
        combat_path = config.combatnorm
        combat_BIN = module_utils.which(combat_path)
        assert combat_BIN, 'cannot find the %s' % combat_path
        command = ['python', combat_BIN, '-f', data_node.identifier, '-o',
                   outfile, '-label', cls_node.identifier]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for combat fails' % outfile
        )
        out_node = bie3.Data(rulebase._SignalFile_Merge, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    return False


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='_SignalFile_Merge')
    filter2 = module_utils.AntecedentFilter(datatype_name='ClassLabelFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def name_outfile(antecedents, user_options):
    data_node, cls_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_combat_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    data_node, cls_node = antecedents
    new_parameters = data_node.data.attributes.copy()
    new_parameters['combat_norm'] = 'yes'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, cls_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
