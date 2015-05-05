#flag_dups_in_bam_folder.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    
    mark_duplicates_path = config.Mark_duplicates
    assert os.path.exists(
        mark_duplicates_path), 'cannot find the %s' % mark_duplicates_path
    for filename in filenames:
        infile = os.path.join(directory, filename)
        outname = os.path.splitext(filename)[-2] + '.bam'
        outname = os.path.join(outfile, outname)
        command = ['java', '-Xmx5g', '-jar', mark_duplicates_path,
                   'I=' + infile, 'O=' + outname, 'METRICS_FILE=metricsFile',
                   'VALIDATION_STRINGENCY=LENIENT', 'REMOVE_DUPLICATES=true']
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()
        if 'error' in error_message[1]:
            raise ValueError(error_message)
        assert module_utils.exists_nz(outname), (
            'the output file %s for flag_dups_in_bam_folder does not exist' %
            outname)
    
    out_node = bie3.Data(rulebase.BamFolder, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'bamFiles_dups' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='BamFolder')
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1)
    return data_node
