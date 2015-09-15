from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from Betsy import module_utils

        fastq_node, sample_node, reference_node = antecedents
        module_utils.safe_mkdir(out_path)

        fastq_path = fastq_node.identifier
        reference_path = reference_node.identifier

        assert os.path.exists(fastq_path)
        assert os.path.exists(reference_path)
        assert os.path.isdir(fastq_path)
        assert os.path.isdir(reference_path)

        bowtie = module_utils.which_assert(config.bowtie)
        reference_genome = module_utils.find_bowtie1_reference(reference_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            sample_node.identifier, fastq_path)
        fastq_files = x

        # Make a list of the jobs to run.
        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            sam_filename = os.path.join(out_path, "%s.sam" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, sam_filename, log_filename
            jobs.append(x)
        
        # Generate bowtie1 commands for each of the files.
        attr2orient = {
            "single" : None,
            "paired_ff" : "ff",
            "paired_fr" : "fr",
            "paired_rf" : "rf",
            }
        x = sample_node.data.attributes["orientation"]
        orientation = attr2orient[x]

        sq = module_utils.shellquote
        commands = []
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = module_utils.make_bowtie1_command(
                reference_genome, sam_filename, pair1, fastq_file2=pair2,
                orientation=orientation, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        module_utils.run_parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            # Make sure sam file created.
            assert module_utils.exists_nz(sam_filename), \
                   "Missing: %s" % sam_filename
            # Make sure there are some alignments.
            x = open(log_filename).read()
            assert x.find("No alignments") < 0, "No alignments"


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #data_node, group_node = antecedents
        #original_file = module_utils.get_inputid(data_node.identifier)
        #filename = 'Samfolder_' + original_file
        #return filename
        return "alignments.bowtie"


## def concatenate_files(input_files, outfile):
##     """Concatenate multiple files into a single one"""
##     import os
##     import gzip
    
##     with open(outfile, 'w') as outfile:
##         for fname in input_files:
##             if fname.endswith('.gz'):
##                 if os.path.exists(fname):
##                     infile = gzip.open(fname)
##                 else:
##                     # see if the unzipped file exists
##                     fname = os.path.splitext(fname)[0]
##                     if os.path.exists(fname):
##                         infile = open(fname)
##                     else:
##                         raise ValueError('cannot find %s' % fname)
##             else:
##                 assert os.path.exists(fname), 'cannot find %s' % fname
##                 infile = open(fname)
##             for line in infile:
##                 outfile.write(line)


## def concatenate_multiple_line(group_dict, foldername):
##     """given group_dict, concatenate multiple line fastq files under foldername
##        output is a new dictionary in format <sample:[R1_sample,R2_sample]>
##                    or <sample:[sample]>"""
##     import os

##     current_dir = os.getcwd()
##     new_group_dict = {}
##     for sample_name, files in group_dict.iteritems():
##         if len(files) == 1:
##             outfile = os.path.join(current_dir, sample_name + '.fastq')
##             inputfiles = [os.path.join(foldername, x) for x in files[0]]
##             concatenate_files(inputfiles, outfile)
##             new_group_dict[sample_name] = [outfile]
##         elif len(files) == 2:
##             outfile_left = os.path.join(current_dir, sample_name + '_R1.fastq')
##             outfile_right = os.path.join(current_dir, sample_name + '_R2.fastq')
##             inputfiles_left = [os.path.join(foldername, x) for x in files[0]]
##             inputfiles_right = [os.path.join(foldername, x) for x in files[1]]
##             concatenate_files(inputfiles_left, outfile_left)
##             concatenate_files(inputfiles_right, outfile_right)
##             new_group_dict[sample_name] = [outfile_left, outfile_right]
##         else:
##             raise ValueError(
##                 'the number of files for each sample is not correct')
##     return new_group_dict


## def preprocess_multiple_sample(folder, group_dict, outfile, ref):
##     """folder:the folder where files are stored,
##        group_dict: A dictionary in format <sample:[[R1_samples],[R2_samples]]>
##                    or <sample:[[samples]]>
##        outfile: output file name,
##        ref: reference species, human or mouse"""
##     import os
##     import subprocess
##     from genomicode import config

##     if not os.path.exists(outfile):
##         os.mkdir(outfile)
##     if ref == 'human':
##         ref_file = config.rna_human
##     elif ref == 'mouse':
##         ref_file = config.rna_mouse
##     else:
##         raise ValueError("we cannot handle %s" % ref)
##     new_group_dict = concatenate_multiple_line(group_dict, folder)
##     for sample in new_group_dict:
##         files = new_group_dict[sample]
##         if len(files) == 1:
##             input_file = os.path.join(folder, files[0])
##             command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
##                        '99999999', '-l', '25', '-p', '8', '-a', '-m', '200',
##                        '-S', ref_file, input_file]
##         elif len(files) == 2:
##             input_file1 = os.path.join(folder, files[0])
##             input_file2 = os.path.join(folder, files[1])
##             command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
##                        '99999999', '-l', '25', '-I', '1', '-X', '1000', '-p',
##                        '8', '-a', '-m', '200', '-S', ref_file, '-1',
##                        input_file1, '-2', input_file2]
##         else:
##             raise ValueError('number files is not correct')
##         outfilename = os.path.join(outfile, sample + '.sam')
##         f = file(outfilename, 'w')
##         try:
##             process = subprocess.Popen(
##                 command, shell=False, stdout=f, stderr=subprocess.PIPE)
##             process.wait()
##             error_message = process.communicate()[1]
##             if 'error' in error_message:
##                 raise ValueError(error_message)
##         finally:
##             f.close()
