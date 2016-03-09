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
        from Betsy import module_utils as mlib

        fastq_node, sample_node, strand_node, reference_node, gene_node = \
                    antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        ref = alignlib.create_reference_genome(reference_node.identifier)
        assert os.path.exists(ref.fasta_file_full)
        gtf_file = gene_node.identifier
        filelib.assert_exists_nz(gtf_file)
        stranded = mlib.read_stranded(strand_node.identifier)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "TopHat %s" % alignlib.get_tophat_version()

        # Get the GTF file, if any.
        #gtf_file = module_utils.get_user_option(
        #    user_options, "tophat_gtf_file", check_file=True)
        transcriptome_fa = mlib.get_user_option(
            user_options, "tophat_transcriptome_fa", check_file=True)
        assert gtf_file or transcriptome_fa, (
            "Either tophat_gtf_file or tophat_transcriptome_fa (preferred) "
            "must be provided.")

        # Make a list of the jobs to run.
        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            tophat_path = os.path.join(out_path, "%s.tophat" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, tophat_path, log_filename
            jobs.append(x)
        
        # Generate tophat commands for each of the files.
        s2ltype = {
            "unstranded" : "fr-unstranded",
            "firststrand" : "fr-firststrand",
            "secondstrand" : "fr-secondstrand",
            }
        assert stranded.stranded in s2ltype, "Unknown stranded: %s" % \
               stranded.stranded
        library_type = s2ltype[stranded.stranded]

        sq = parallel.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, tophat_path, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_tophat_command(
                ref.fasta_file_full, tophat_path, pair1, fastq_file2=pair2,
                gtf_file=gtf_file, transcriptome_fa=transcriptome_fa,
                library_type=library_type, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[3] for x in jobs]  # out_path
        x = [os.path.join(x, "accepted_hits.bam") for x in x]
        bam_filenames = x
        filelib.assert_exists_nz_many(bam_filenames)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "alignments.tophat"
