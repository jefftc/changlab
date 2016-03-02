from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils

        fastq_node, sample_node, reference_node = antecedents
        fastq_path = fastq_node.identifier
        ref = alignlib.create_reference_genome(reference_node.identifier)
        filelib.safe_mkdir(out_path)

        reference_path = reference_node.identifier

        assert os.path.exists(fastq_path)
        assert os.path.exists(ref.fasta_file_full)
        assert os.path.isdir(fastq_path)

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
        
        # Generate bowtie2 commands for each of the files.
        attr2orient = {
            "single" : None,
            "paired" : None,
            "paired_ff" : "ff",
            "paired_fr" : "fr",
            "paired_rf" : "rf",
            }
        x = sample_node.data.attributes["orientation"]
        orientation = attr2orient[x]

        sq = parallel.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_bowtie2_command(
                ref.fasta_file_full, pair1, fastq_file2=pair2,
                orientation=orientation, sam_file=sam_filename, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            assert filelib.exists_nz(sam_filename), \
                   "Missing: %s" % sam_filename


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #data_node, group_node = antecedents
        #original_file = module_utils.get_inputid(data_node.identifier)
        #filename = 'Samfolder_' + original_file
        #return filename
        return "alignments.bowtie2"
