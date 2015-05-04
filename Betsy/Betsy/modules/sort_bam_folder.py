#sort_bam_folder.py
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
    
    samtools_BIN = config.samtools
    assert os.path.exists(samtools_BIN), 'cannot find the %s' % samtools_BIN
    for filename in filenames:
        infile = os.path.join(directory, filename)
        outname = os.path.splitext(filename)[-2] + '_sorted.bam'
        outname = os.path.join(outfile, outname)
        command = [samtools_BIN, 'sort', infile, outname]

        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()
        if 'error' in error_message[1]:
            raise ValueError(error_message)
        assert module_utils.exists_nz(outname), (
            'the output file %s for sort_bam_folder does not exist' % outname
        )
    
    out_node = bie3.Data(rulebase.BamFolder, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'bamFiles_sorted' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)

    return data_node
