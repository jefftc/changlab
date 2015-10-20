from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import alignlib
        from genomicode import shell
        from Betsy import module_utils
        
        fastq_node, sample_node, reference_node = antecedents
        fastq_path = fastq_node.identifier
        assert os.path.exists(fastq_path)
        assert os.path.isdir(fastq_path)
        ref = alignlib.create_reference_genome(reference_node.identifier)
        filelib.safe_mkdir(out_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            sample_node.identifier, fastq_path)
        fastq_files = x
        assert fastq_files, "I could not find any FASTQ files."

        # Make a list of the jobs to run.
        jobs = []  # list of sample, pair1, pair2, log_filename
        for x in fastq_files:
            sample, pair1, pair2 = x
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, log_filename
            jobs.append(x)
        
        sq = shell.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_rsem_command(
                ref.fasta_file_full, sample, pair1, fastq_file2=pair2,
                num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        # Need to run in out_path.  Otherwise, files will be everywhere.
        shell.parallel(commands, max_procs=num_cores, path=out_path)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1, pair2, log_filename = x
            # Make sure sam file created.
            assert filelib.exists_nz(log_filename), \
                   "Missing: %s" % log_filename
            # Make sure there are some alignments.
            filename = os.path.join(out_path, "%s.genes.results" % sample)
            assert filelib.exists_nz(filename)

        
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #data_node, group_node = antecedents
        #original_file = module_utils.get_inputid(data_node.identifier)
        #filename = original_file + '.tdf'
        #return filename
        return "rsem"

## def preprocess_single_sample(folder, sample, files, out_file, ref, num_cores):
##     import os
##     import subprocess
##     from Betsy import module_utils
##     from genomicode import config
##     if ref == 'human':
##         ref_file = config.rna_human
##     elif ref == 'mouse':
##         ref_file = config.rna_mouse
##     else:
##         raise ValueError("we cannot handle %s" % ref)
    
##     bamfiles = os.listdir(folder)
##     if sample + '.bam' in bamfiles:
##         input_file = os.path.join(folder, sample + '.bam')
##         command = ['rsem-calculate-expression', '--bam', input_file, ref_file,
##                    '--no-bam-output', '-p', str(num_cores), sample]
##         if len(files) == 2:
##             command.append('--paired-end')
##     else:
##         if len(files) == 1:
##             input_file = os.path.join(folder, files[0][0])
##             command = ['rsem-calculate-expression', '--bam', input_file,
##                        ref_file, '--no-bam-output', '-p', str(num_cores),
##                        sample]
##         elif len(files) == 2:
##             input_file1 = os.path.join(folder, files[0][0])
##             input_file2 = os.path.join(folder, files[1][0])
##             command = ['rsem-calculate-expression', '--bam', '--paired-end',
##                        input_file1, input_file2, ref_file, '--no-bam-output',
##                        '-p', str(num_cores), sample]
##         else:
##             raise ValueError('number files is not correct')
    
##     process = subprocess.Popen(command,
##                                shell=False,
##                                stdout=subprocess.PIPE,
##                                stderr=subprocess.PIPE)
##     #process.wait()
##     error_message = process.communicate()[1]
##     if 'error' in error_message:
##         raise ValueError(error_message)
    
##     slice_BIN = config.slice_matrix
##     command = ['python', slice_BIN, sample + '.genes.results',
##                '--select_col_ids', 'transcript_id,gene_id,TPM',
##                '--replace_col_ids', 'TPM,' + sample]
##     f = file(out_file, 'w')
##     try:
##         process = subprocess.Popen(command,
##                                    shell=False,
##                                    stdout=f,
##                                    stderr=subprocess.PIPE)
##         process.wait()
##     finally:
##         f.close()
    
##     error_message = process.communicate()[1]
##     if 'error' in error_message:
##         raise ValueError(error_message)
    
##     assert module_utils.exists_nz(out_file), (
##         'the output file %s does not exist' % out_file
##     )


## def preprocess_multiple_sample(folder, group_dict, outfile, ref, num_cores):
##     import os
##     import tempfile
##     from Betsy import module_utils
    
##     #filenames = os.listdir(folder)
##     file_list = []
##     for sample in group_dict:
##         temp_file = tempfile.mkstemp()[1]
##         preprocess_single_sample(
##             folder, sample, group_dict[sample], temp_file,
##             ref, num_cores)
##         file_list.append(temp_file)
    
##     result_file = file_list[0]
##     tmp_list = file_list[:]
##     try:
##         for filename in file_list[1:]:
##             tmp_result = tempfile.mkstemp()[1]
##             f = file(tmp_result, 'a+')
##             try:
##                 module_utils.merge_two_files(result_file, filename, f)
##             finally:
##                 f.close()
##                 tmp_list.append(tmp_result)
##             result_file = tmp_result
##         os.rename(result_file, outfile)
##     finally:
##         for filename in tmp_list:
##             if os.path.exists(filename):
##                 os.remove(filename)
    
