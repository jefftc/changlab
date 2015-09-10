from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from Betsy import module_utils
        from genomicode import config

        fastq_node, group_node, index_node = antecedents
        module_utils.safe_mkdir(out_path)

        fastq_path = fastq_node.identifier
        index_path = index_node.identifier

        assert os.path.exists(fastq_path)
        assert os.path.exists(index_path)
        assert os.path.isdir(fastq_path)
        assert os.path.isdir(index_path)

        # bwa aln -t <num_cores> <reference.fa> <input.fastq> > <output.sai>
        bwa = module_utils.which_assert(config.bwa)

        # Find the indexed reference genome at:
        # <index_path>/<assembly>.fa
        x = os.listdir(index_path)
        x = [x for x in x if x.lower().endswith(".fa")]
        assert len(x) == 1, "Cannot find bwa index."
        x = x[0]
        reference_fa = os.path.join(index_path, x)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            group_node.identifier, fastq_path)
        grouped_fastq_files = x

        # Make sure no duplicate samples.
        x1 = [x[0] for x in grouped_fastq_files]
        x2 = {}.fromkeys(x1).keys()
        assert len(x1) == len(x2), "dup sample"

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
            x = x + ".sai"
            out_file = x
            out_filename = os.path.join(out_path, out_file)
            x = in_filename, out_filename
            jobs.append(x)

        # Calculate the number of cores per job.
        nc = 1
        nc = max(1, num_cores/len(jobs))

        # Make the bwa commands.
        sq = module_utils.shellquote
        commands = []
        for x in jobs:
            fastq_filename, sai_filename = x

            x = [
                bwa,
                "aln",
                "-t", nc,
                sq(reference_fa),
                sq(fastq_filename),
                ">",
                sq(sai_filename),
                ]
            x = " ".join(map(str, x))
            commands.append(x)

        module_utils.run_parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            in_filename, out_filename = x
            assert module_utils.exists_nz(out_filename), \
                   "Missing: %s" % out_filename

        
    def name_outfile(self, antecedents, user_options):
        return "alignments.bwa.sai"

