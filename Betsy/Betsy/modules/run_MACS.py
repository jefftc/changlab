from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from Betsy import module_utils

        bam_node, group_node, index_node = antecedents
        module_utils.safe_mkdir(out_path)

        bam_path = fastq_node.identifier
        assert os.path.exists(bam_path)
        assert os.path.isdir(bam_path)

        # Find the BAM files.
        raise NotImplementedError


        # Make a list of all the jobs to do.
        jobs = []   # list of (sample, pair1, pair2, sam_filename)
        for x in grouped_fastq_files:
            sample, pair1, pair2 = x
            sam_filename = os.path.join(out_path, "%s.sam" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, sam_filename, log_filename
            jobs.append(x)

        # Calculate the number of cores per job.
        nc = max(1, num_cores/len(jobs))
        
        # Make the bwa commands.
        commands = []
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            x = module_utils.make_bwa_mem_command(
                reference_fa, sam_filename, log_filename, pair1,
                fastq_file2=pair2, num_threads=nc)
            commands.append(x)

        module_utils.run_parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            assert module_utils.exists_nz(sam_filename), \
                   "Missing: %s" % sam_filename
            # TODO: Should also check log file.

        
    def name_outfile(self, antecedents, user_options):
        return "MACS"


