from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import shell
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

        # Get the GTF file, if any.
        gtf_file = module_utils.get_user_option(
            user_options, "tophat_gtf_file", check_file=True)
        transcriptome_fa = module_utils.get_user_option(
            user_options, "tophat_transcriptome_fa", check_file=True)
        assert gtf_file or transcriptome_fa, \
               ("Either tophat_gtf_file or tophat_transcriptome_index must be "
                "provided.")

        # Make a list of the jobs to run.
        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            tophat_path = os.path.join(out_path, "%s.tophat" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, tophat_path, log_filename
            jobs.append(x)
        
        # Generate tophat commands for each of the files.
        attr2orient = {
            "single" : None,
            "paired" : None,
            "paired_ff" : "ff",
            "paired_fr" : "fr",
            "paired_rf" : "rf",
            }
        x = sample_node.data.attributes["orientation"]
        orientation = attr2orient[x]

        sq = shell.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, tophat_path, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_tophat_command(
                ref.fasta_file_full, tophat_path, pair1, fastq_file2=pair2,
                gtf_file=gtf_file, transcriptome_fa=transcriptome_fa,
                orientation=orientation, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[3] for x in jobs]  # out_path
        x = [os.path.join(x, "accepted_hits.bam") for x in x]
        bam_filenames = x
        filelib.assert_exists_nz_many(bam_filenames)


    def name_outfile(self, antecedents, user_options):
        return "alignments.tophat"
