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

        fastq_node, sample_node, ref_node = antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # Do a quick check to make sure the reference is correct.
        # Otherwise, error may be hard to disgnose.
        alignlib.assert_is_STAR_reference(ref.path)

        metadata = {}
        metadata["tool"] = "STAR %s" % alignlib.get_STAR_version()

        # Figure out the strandedness.
        is_stranded = False

        # STAR --runThreadN 40 --genomeDir test05 \
        #   --readFilesIn test.fastq/test03_R1_001.fastq \
        #   test.fastq/test03_R2_001.fastq --outFileNamePrefix test06.
        # If unstranded, add --outSAMstrandField intronMotif
        
        # Make a list of the jobs to run.
        jobs = []  # list of filelib.GenericObject objects
        for x in fastq_files:
            sample, pair1, pair2 = x
            out_prefix = "%s." % sample
            bam_filename = os.path.join(
                out_path, "%sAligned.out.bam" % out_prefix)
            log_filename = os.path.join(out_path, "%s.log" % sample)

            x = filelib.GenericObject(
                sample=sample, pair1=pair1, pair2=pair2,
                out_prefix=out_prefix,
                bam_filename=bam_filename,
                log_filename=log_filename,
                )
            jobs.append(x)

        # Run pass 1.
        commands = []
        for j in jobs:
            x = os.path.join(out_path, j.out_prefix)
            cmd = alignlib.make_STAR_command(
                ref.path, x, num_cores, is_stranded, j.pair1, j.pair2,
                j.log_filename)
            # For debugging.  If this file already exists, skip it.
            if not filelib.exists_nz(j.bam_filename):
                parallel.sshell(cmd, path=out_path)
            filelib.assert_exists_nz(j.bam_filename)
            commands.append(cmd)
        
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        
        return metadata

    def name_outfile(self, antecedents, user_options):
        return "alignments.star"
