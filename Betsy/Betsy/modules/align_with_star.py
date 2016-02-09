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
        from genomicode import config
        from Betsy import module_utils

        fastq_node, sample_node, reference_node = antecedents
        fastq_path = fastq_node.identifier
        reference_path = reference_node.identifier
        assert os.path.exists(fastq_path)
        assert os.path.isdir(fastq_path)
        assert os.path.exists(reference_path)
        assert os.path.isdir(reference_path)
        filelib.safe_mkdir(out_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            sample_node.identifier, fastq_path)
        fastq_files = x

        # Figure out the orientation.
        x = sample_node.data.attributes["orientation"]
        assert x == "single" or x.startswith("paired")
        is_paired = x.startswith("paired")
        is_stranded = False
        if is_paired:
            is_stranded = (x != "paired")

        # STAR --runThreadN 40 --genomeDir test05 \
        #   --readFilesIn test.fastq/test03_R1_001.fastq \
        #   test.fastq/test03_R2_001.fastq --outFileNamePrefix test06.
        # If unstranded, add --outSAMstrandField intronMotif
        
        STAR = filelib.which_assert(config.STAR)

        # Make a list of the jobs to run.
        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            out_prefix = "%s." % sample
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, out_prefix, log_filename
            jobs.append(x)

        # TODO: Play around with runThreadN parameter.

        # Make the commands.
        sq = shell.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, out_prefix, log_filename = x
            
            x = [
                sq(STAR),
                "--genomeDir", sq(reference_path),
                "--outFileNamePrefix", out_prefix,
                ]
            if not is_stranded:
                x += ["--outSAMstrandField", "intronMotif"]
            x += ["--readFilesIn", sq(pair1)]
            if pair2:
                x += [sq(pair2)]
            x = " ".join(map(str, x))
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        # STAR takes 27 Gb per process.  Make sure we don't use up
        # more memory than is available on the machine.
        max_procs = module_utils.calc_max_procs_from_ram(30)
        nc = min(num_cores, max_procs)
        shell.parallel(commands, max_procs=nc, path=out_path)

        # Make sure the analysis completed successfully.
        x = [x[-2] for x in jobs]  # out_prefix
        x = ["%sAligned.out.sam" % x for x in x]
        x = [os.path.join(out_path, x) for x in x]
        filelib.assert_exists_nz_many(x)


    def name_outfile(self, antecedents, user_options):
        return "alignments.star"
