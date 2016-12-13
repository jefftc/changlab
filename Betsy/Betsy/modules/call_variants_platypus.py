from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib
        import call_variants_GATK

        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # Figure out whether the user wants SNPs or INDELs.
        #assert "vartype" in out_attributes
        #vartype = out_attributes["vartype"]
        #assert vartype in ["all", "snp", "indel"]

        # Platypus generates an error if there are spaces in the BAM
        # filename.  Symlink the file to a local directory to make
        # sure there are no spaces.
        bam_path = "bam"

        jobs = []  # list of filelib.GenericObject
        for bam_filename in bam_filenames:
            p, f = os.path.split(bam_filename)
            sample, ext = os.path.splitext(f)
            bai_filename = "%s.bai" % bam_filename
            filelib.assert_exists_nz(bai_filename)
            x = sample.replace(" ", "_")
            local_bam = os.path.join(bam_path, "%s.bam" % x)
            local_bai = os.path.join(bam_path, "%s.bam.bai" % x)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            err_filename = os.path.join(out_path, "%s.err" % sample)
            # Unfiltered file.
            #raw_filename = os.path.join(out_path, "%s.raw" % sample)
            # Final VCF file.
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = filelib.GenericObject(
                bam_filename=bam_filename, bai_filename=bai_filename,
                local_bam=local_bam, local_bai=local_bai,
                log_filename=log_filename, err_filename=err_filename,
                out_filename=out_filename)
            jobs.append(x)

        filelib.safe_mkdir(bam_path)
        for j in jobs:
            assert " " not in j.local_bam
            filelib.assert_exists_nz(j.bam_filename)
            filelib.assert_exists_nz(j.bai_filename)
            if not os.path.exists(j.local_bam):
                os.symlink(j.bam_filename, j.local_bam)
            if not os.path.exists(j.local_bai):
                os.symlink(j.bai_filename, j.local_bai)

        # TODO: Keep better track of the metadata.
        buffer_size = 100000
        max_reads = 5E6
        # Running into errors sometimes, so increase these numbers.
        #   WARNING - Too many reads (5000000) in region
        #   1:500000-600000. Quitting now. Either reduce --bufferSize or
        #   increase --maxReads.
        buffer_size = buffer_size * 10
        max_reads = max_reads * 10
        
        # Make a list of commands.
        commands = []
        for j in jobs:
            #nc = max(1, num_cores/len(jobs))
            x = alignlib.make_platypus_command(
                bam_file=j.local_bam,
                ref_file=ref.fasta_file_full,
                log_file=j.log_filename,
                out_file=j.out_filename,
                buffer_size=buffer_size, max_reads=max_reads)
            x = "%s >& %s" % (x, j.err_filename)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.  If not, try
        # to diagnose.
        for j in jobs:
            if filelib.exists_nz(j.out_filename):
                continue
            for line in open(j.err_filename):
                if line.find("WARNING - Too many reads") >= 0:
                    print line,
        x = [j.out_filename for j in jobs]
        filelib.assert_exists_nz_many(x)

        # Filter each of the VCF files.
        #for j in jobs:
        #    call_variants_GATK.filter_by_vartype(
        #        vartype, j.raw_filename, j.out_filename)
        #metadata["filter"] = vartype

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "platypus.vcf"
