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

        bam_node, nc_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        metadata["tool"] = "MuSE %s" % alignlib.get_muse_version()

        wgs_or_wes = mlib.get_user_option(
            user_options, "wgs_or_wes", not_empty=True,
            allowed_values=["wgs", "wes"])
        dbsnp_file = mlib.get_user_option(
            user_options, "muse_dbsnp_vcf", not_empty=True, check_file=True)

        # Make sure dbsnp_file is compressed and indexed.
        assert dbsnp_file.endswith(".vcf.gz"), \
               "muse_dbsnp_vcf must be bgzip compressed."
        x = "%s.tbi" % dbsnp_file
        assert filelib.exists_nz(x), "muse_dbsnp_vcf must be tabix indexed."

        # sample -> bam filename
        sample2bamfile = mlib.root2filename(bam_filenames)
        # Make sure files exist for all the samples.
        mlib.assert_normal_cancer_samples(nc_match, sample2bamfile)

        # list of (normal_sample, cancer_sample, normal_bamfile, tumor_bamfile,
        #   muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile,
        #   logfile1, logfile2)
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            muse_call_stem = opj(out_path, "%s.call" % cancer_sample)
            muse_call_file = "%s.MuSE.txt" % muse_call_stem
            raw_vcf_outfile = opj(out_path, "%s.vcf.raw" % cancer_sample)
            vcf_outfile = opj(out_path, "%s.vcf" % cancer_sample)
            log_outfile1 = opj(out_path, "%s.call.log" % cancer_sample)
            log_outfile2 = opj(out_path, "%s.sump.log" % cancer_sample)
            x = normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile, \
                log_outfile1, log_outfile2
            jobs.append(x)

        # Generate the commands.
        # MuSE call -O test11 -f genomes/Broad.hg19/Homo_sapiens_assembly19.fa\
        #   bam04/196B-MG.bam bam04/PIM001_G.bam
        # MuSE sump -I test11.MuSE.txt -E -O test12.vcf \
        #   -D MuSE/dbsnp_132_b37.leftAligned.vcf.gz

        MuSE = mlib.findbin("muse")
        
        sq = mlib.sq
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile, \
                log_outfile1, log_outfile2 = x

            x = [
                sq(MuSE),
                "call",
                "-O", muse_call_stem,
                "-f", sq(ref.fasta_file_full),
                cancer_bamfile,
                normal_bamfile,
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, log_outfile1)
            commands.append(x)
        assert len(commands) == len(jobs)
        # Not sure about RAM.
        nc = mlib.calc_max_procs_from_ram(10, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["num_cores"] = nc
        metadata["commands"] = commands

        # Make sure the log files have no errors.  The files should be
        # empty.
        log_files = [x[8] for x in jobs]
        filelib.assert_exists_z_many(log_files)

        # Make sure the call files are created and not empty.
        call_files = [x[5] for x in jobs]
        filelib.assert_exists_nz_many(call_files)

        # Run the "sump" step.
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile, \
                log_outfile1, log_outfile2 = x

            x = [
                sq(MuSE),
                "sump",
                "-I", sq(muse_call_file),
                ]
            assert wgs_or_wes in ["wgs", "wes"]
            if wgs_or_wes == "wgs":
                x += ["-G"]
            else:
                x += ["-E"]
            x += [
                "-O", sq(raw_vcf_outfile),
                "-D", sq(dbsnp_file),
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, log_outfile2)
            commands.append(x)
        assert len(commands) == len(jobs)
        # Not sure about RAM.
        nc = mlib.calc_max_procs_from_ram(10, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["commands"] = metadata["commands"] + commands

        # Make sure the log files have no errors.  The files should be
        # empty.
        log_files = [x[9] for x in jobs]
        filelib.assert_exists_z_many(log_files)

        # Make sure the raw files are created and not empty.
        vcf_files = [x[6] for x in jobs]
        filelib.assert_exists_nz_many(vcf_files)


        # Fix the files.
        commands = []  # Should be python commands.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile, \
                log_outfile1, log_outfile2 = x
            args = normal_sample, cancer_sample, raw_vcf_outfile, vcf_outfile
            x = alignlib.clean_muse_vcf, args, {}
            commands.append(x)
        parallel.pyfun(commands, num_procs=num_cores)
        
        # Delete the log_outfiles if empty.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                muse_call_stem, muse_call_file, raw_vcf_outfile, vcf_outfile, \
                log_outfile1, log_outfile2 = x
            if os.path.exists(log_outfile1):
                os.unlink(log_outfile1)
            if os.path.exists(log_outfile2):
                os.unlink(log_outfile2)

        # Make sure output VCF files exist.
        x = [x[7] for x in jobs]
        filelib.assert_exists_many(x)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "muse.vcf"
