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

        fastq_node, sample_node, strand_node, ref_node = antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        stranded = mlib.read_stranded(strand_node.identifier)
        filelib.safe_mkdir(out_path)

        # Do a quick check to make sure the reference is correct.
        # Otherwise, error may be hard to disgnose.
        alignlib.assert_is_STAR_reference(ref.path)
        

        metadata = {}
        metadata["tool"] = "STAR %s" % alignlib.get_STAR_version()

        x = mlib.get_user_option(
            user_options, "two_pass", allowed_values=["no", "yes"])
        two_pass = (x == "yes")

        # Figure out the strandedness.
        is_stranded = stranded.stranded != "unstranded"

        # STAR --runThreadN 40 --genomeDir test05 \
        #   --readFilesIn test.fastq/test03_R1_001.fastq \
        #   test.fastq/test03_R2_001.fastq --outFileNamePrefix test06.
        # If unstranded, add --outSAMstrandField intronMotif
        
        # Make a list of the jobs to run.
        jobs = []  # list of filelib.GenericObject objects
        for x in fastq_files:
            sample, pair1, pair2 = x
            pass1_out_prefix = "p1.%s." % sample
            pass2_out_prefix = "%s." % sample
            pass1_bam_filename = os.path.join(
                out_path, "%sAligned.out.bam" % pass1_out_prefix)
            pass2_bam_filename = os.path.join(
                out_path, "%sAligned.out.bam" % pass2_out_prefix)
            sjdb_filename = os.path.join(out_path, "p1.%s.SJ.out.tab" % sample)
            log1_filename = os.path.join(out_path, "p1.%s.log" % sample)
            log2_filename = os.path.join(out_path, "%s.log" % sample)

            x = filelib.GenericObject(
                sample=sample, pair1=pair1, pair2=pair2,
                pass1_out_prefix=pass1_out_prefix, 
                pass2_out_prefix=pass2_out_prefix,
                pass1_bam_filename=pass1_bam_filename,
                pass2_bam_filename=pass2_bam_filename,
                sjdb_filename=sjdb_filename,
                log1_filename=log1_filename,
                log2_filename=log2_filename,
                )
            jobs.append(x)

        # Run pass 1.
        commands = []
        for j in jobs:
            x = os.path.join(out_path, j.pass1_out_prefix)
            cmd = alignlib.make_STAR_command(
                ref.path, x, num_cores, is_stranded, j.pair1, j.pair2,
                j.log1_filename)
            # For debugging.  If this file already exists, skip it.
            if not filelib.exists_nz(j.pass1_bam_filename):
                parallel.sshell(cmd, path=out_path)
            filelib.assert_exists_nz(j.pass1_bam_filename)
            commands.append(cmd)

        if two_pass:
            # Make a new index with the splice junction information.
            sj_index = os.path.join(out_path, "genome.2pass")
            x = [x.sjdb_filename for x in jobs]
            filelib.assert_exists_nz_many(x)
            x = alignlib.make_STAR_index_command(
                ref.fasta_file_full, sj_index, sjdb_files=x,
                num_cores=num_cores)
            x = "%s >& genome.2pass.log" % x
            commands.append(x)

            # For debugging.  If this file already exists, skip it.
            if filelib.exists_nz("genome.2pass.log"):
                parallel.sshell(x, path=out_path)
            alignlib.assert_is_STAR_reference(sj_index)

        # Run pass 2.
        for j in jobs:
            # For debugging.  If this file already exists, skip it.
            if os.path.exists(j.pass2_bam_filename):
                continue
            if two_pass:
                x = os.path.join(out_path, j.pass2_out_prefix)
                cmd = alignlib.make_STAR_command(
                    sj_index, x, num_cores, is_stranded, j.pair1, j.pair2,
                    j.log2_filename)
                parallel.sshell(cmd, path=out_path)
                commands.append(cmd)
            else:
                # link pass1_bam_filename to pass2_bam_filename
                os.symlink(j.pass1_bam_filename, j.pass2_bam_filename)
                continue
            filelib.assert_exists_nz(j.pass2_bam_filename)
        
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        
        # STAR takes 28 Gb per process.  Make sure we don't use up
        # more memory than is available on the machine.
        # Defaults:
        # --limitGenomeGenerateRAM   31000000000
        # --outFilterMismatchNmax    10             Num mismatches.
        #nc = mlib.calc_max_procs_from_ram(50, buffer=100, upper_max=num_cores)
        #metadata["num_cores"] = nc
        #parallel.pshell(commands, max_procs=nc, path=out_path)
        
        # Make sure the analysis completed successfully.
        #x = [x[-2] for x in jobs]  # sam_filename
        #filelib.assert_exists_nz_many(x)
        return metadata

    def name_outfile(self, antecedents, user_options):
        return "alignments.star"
            



