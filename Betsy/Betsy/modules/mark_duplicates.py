#mark_duplcates.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    mark_duplicates_path = config.Mark_duplicates
    assert os.path.exists(
        mark_duplicates_path), 'cannot find the %s' % mark_duplicates_path
    command = ['java', '-Xmx5g', '-jar', mark_duplicates_path,
               'I=' + in_data.identifier, 'O=' + outfile,
               'METRICS_FILE=metricsFile', 'VALIDATION_STRINGENCY=LENIENT',
               'REMOVE_DUPLICATES=true']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    #error_message = process.communicate()[1]
    #if 'error' in error_message:
    #    raise ValueError(error_message)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for mark_duplcates does not exist' % outfile
    )
    out_node = bie3.Data(rulebase.BamFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'marked_duplicates_' + original_file + '.bam'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='BamFile')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node
