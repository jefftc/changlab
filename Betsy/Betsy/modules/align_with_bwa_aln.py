from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils

        fastq_node, group_node, reference_node = antecedents
        fastq_path = fastq_node.identifier
        assert os.path.exists(fastq_path)
        assert os.path.isdir(fastq_path)
        ref = alignlib.create_reference_genome(reference_node.identifier)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "bwa %s" % alignlib.get_bwa_version()

        #reference_fa = module_utils.find_bwa_reference(index_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            group_node.identifier, fastq_path)
        grouped_fastq_files = x

        # Make sure no duplicate samples.
        x1 = [x[0] for x in grouped_fastq_files]
        x2 = {}.fromkeys(x1).keys()
        assert len(x1) == len(x2), "dup sample"

        # Make a list of all FASTQ files to align.
        fastq_files = []
        for x in grouped_fastq_files:
            sample, pair1, pair2 = x
            assert pair1
            fastq_files.append(pair1)
            if pair2:
                fastq_files.append(pair2)

        # Make a list of all the jobs to do.
        jobs = []   # list of (fastq_filename, sai_filename)
        for in_filename in fastq_files:
            in_path, in_file = os.path.split(in_filename)
            x = in_file
            if x.lower().endswith(".fq"):
                x = x[:-3]
            elif x.lower().endswith(".fastq"):
                x = x[:-6]
            sai_filename = os.path.join(out_path, "%s.sai" % x)
            log_filename = os.path.join(out_path, "%s.log" % x)
            x = in_filename, sai_filename, log_filename
            jobs.append(x)

        # Calculate the number of cores per job.
        nc = max(1, num_cores/len(jobs))
        metadata["num cores"] = nc

        # Make the bwa commands.
        commands = []
        for x in jobs:
            fastq_filename, sai_filename, log_filename = x
            x = alignlib.make_bwa_aln_command(
                ref.fasta_file_full, fastq_filename, sai_filename,
                log_filename, num_threads=nc)
            commands.append(x)
        metadata["commands"] = commands
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            in_filename, sai_filename, log_filename = x
            assert filelib.exists_nz(sai_filename), \
                   "Missing: %s" % sai_filename
        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "alignments.bwa.sai"
