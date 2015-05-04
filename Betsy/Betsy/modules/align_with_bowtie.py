#align_with_bowtie.py
import os
from Betsy import module_utils, bie3, rulebase
import subprocess
from genomicode import config
import gzip


def concatenate_files(input_files, outfile):
    """Concatenate multiple files into a single one"""
    with open(outfile, 'w') as outfile:
        for fname in input_files:
            if fname.endswith('.gz'):
                if os.path.exists(fname):
                    infile = gzip.open(fname)
                else:
                    # see if the unzipped file exists
                    fname = os.path.splitext(fname)[0]
                    if os.path.exists(fname):
                        infile = open(fname)
                    else:
                        raise ValueError('cannot find %s' % fname)
            else:
                assert os.path.exists(fname), 'cannot find %s' % fname
                infile = open(fname)
            for line in infile:
                outfile.write(line)


def concatenate_multiple_line(group_dict, foldername):
    """given group_dict, concatenate multiple line fastq files under foldername
       output is a new dictionary in format <sample:[R1_sample,R2_sample]>
                   or <sample:[sample]>"""
    current_dir = os.getcwd()
    new_group_dict = {}
    for sample_name, files in group_dict.iteritems():
        if len(files) == 1:
            outfile = os.path.join(current_dir, sample_name + '.fastq')
            inputfiles = [os.path.join(foldername, x) for x in files[0]]
            concatenate_files(inputfiles, outfile)
            new_group_dict[sample_name] = [outfile]
        elif len(files) == 2:
            outfile_left = os.path.join(current_dir, sample_name + '_R1.fastq')
            outfile_right = os.path.join(current_dir, sample_name + '_R2.fastq')
            inputfiles_left = [os.path.join(foldername, x) for x in files[0]]
            inputfiles_right = [os.path.join(foldername, x) for x in files[1]]
            concatenate_files(inputfiles_left, outfile_left)
            concatenate_files(inputfiles_right, outfile_right)
            new_group_dict[sample_name] = [outfile_left, outfile_right]
        else:
            raise ValueError(
                'the number of files for each sample is not correct')
    return new_group_dict


def preprocess_multiple_sample(folder, group_dict, outfile, ref):
    """folder:the folder where files are stored,
       group_dict: A dictionary in format <sample:[[R1_samples],[R2_samples]]>
                   or <sample:[[samples]]>
       outfile: output file name,
       ref: reference species, human or mouse"""
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    if ref == 'human':
        ref_file = config.rna_human
    elif ref == 'mouse':
        ref_file = config.rna_mouse
    else:
        raise ValueError("we cannot handle %s" % ref)
    new_group_dict = concatenate_multiple_line(group_dict, folder)
    for sample in new_group_dict:
        files = new_group_dict[sample]
        if len(files) == 1:
            input_file = os.path.join(folder, files[0])
            command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
                       '99999999', '-l', '25', '-p', '8', '-a', '-m', '200',
                       '-S', ref_file, input_file]
        elif len(files) == 2:
            input_file1 = os.path.join(folder, files[0])
            input_file2 = os.path.join(folder, files[1])
            command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
                       '99999999', '-l', '25', '-I', '1', '-X', '1000', '-p',
                       '8', '-a', '-m', '200', '-S', ref_file, '-1',
                       input_file1, '-2', input_file2]
        else:
            raise ValueError('number files is not correct')
        outfilename = os.path.join(outfile, sample + '.sam')
        f = file(outfilename, 'w')
        try:
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=f,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()[1]
            if 'error' in error_message:
                raise ValueError(error_message)
        finally:
            f.close()


def run(network, antecedents, out_attributes, user_options, num_cores):
    data_node, group_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    group_dict = module_utils.process_group_info(group_node.identifier)
    preprocess_multiple_sample(data_node.identifier, group_dict, outfile,
                               out_attributes['ref'])
    assert module_utils.exists_nz(outfile), (
        'the output file %s for align_with_bowtie does not exist' % outfile
    )
    out_node = bie3.Data(rulebase.SamFolder, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, group_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    data_node, group_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'Samfolder_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes,
                                            datatype='FastqFolder')
    group_node = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='SampleGroupFile')

    return data_node, group_node
