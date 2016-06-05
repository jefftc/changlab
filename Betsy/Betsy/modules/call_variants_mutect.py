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

        bam_node, nc_node, ref_node, interval_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.assert_exists_nz(interval_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        # TODO: Figure out MuTect version.

        # Make sure intervals file ends with:
        # .bed, .list, .picard, .interval_list, or .intervals
        x, x, ext = mlib.splitpath(interval_node.identifier)
        assert ext in [
            ".bed", ".list", ".picard", ".interval_list", ".intervals"]

        cosmic_file = mlib.get_user_option(
            user_options, "mutect_cosmic_vcf", not_empty=True, check_file=True)
        dbsnp_file = mlib.get_user_option(
            user_options, "mutect_dbsnp_vcf", not_empty=True, check_file=True)

        # sample -> bam filename
        sample2bamfile = mlib.root2filename(bam_filenames)
        # Make sure files exist for all the samples.
        mlib.assert_normal_cancer_samples(nc_match, sample2bamfile)

        # list of (cancer_sample, normal_bamfile, tumor_bamfile, call_outfile,
        #    coverage_outfile, vcf_outfile, logfile)
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            call_outfile = opj(out_path, "%s.call_stats.out" % sample)
            cov_outfile = opj(out_path, "%s.coverage.wig.txt" % sample)
            raw_vcf_outfile = opj(out_path, "%s.raw.vcf" % sample)
            vcf_outfile = opj(out_path, "%s.vcf" % sample)
            log_outfile = opj(out_path, "%s.log" % sample)
            x = normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, raw_vcf_outfile, vcf_outfile, \
                log_outfile
            jobs.append(x)

        # java -Xmx2g -jar muTect.jar
        #   --analysis_type MuTect
        #   --reference_sequence <reference>
        #   --cosmic <cosmic.vcf>
        #   --dbsnp <dbsnp.vcf>
        #   --intervals <intervals_to_process>
        #   --input_file:normal <normal.bam>
        #   --input_file:tumor <tumor.bam>
        #   --out <call_stats.out>
        #   --coverage_file <coverage.wig.txt>

        # Generate the commands.
        sq = mlib.sq
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, raw_vcf_outfile, vcf_outfile, \
                log_outfile = x

            UNHASHABLE = [
                ("input_file:normal", sq(normal_bamfile)),
                ("input_file:tumor", sq(cancer_bamfile)),
                ]
            x = alignlib.make_MuTect_command(
                analysis_type="MuTect",
                reference_sequence=sq(ref.fasta_file_full),
                cosmic=sq(cosmic_file),
                dbsnp=sq(dbsnp_file),
                intervals=sq(interval_node.identifier),
                out=sq(call_outfile),
                coverage_file=sq(cov_outfile),
                vcf=sq(raw_vcf_outfile),
                _UNHASHABLE=UNHASHABLE,
                )
            x = "%s >& %s" % (x, log_outfile)
            commands.append(x)
        assert len(commands) == len(jobs)
        nc = mlib.calc_max_procs_from_ram(15, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["num_cores"] = nc
        metadata["commands"] = commands

        # Make sure log files have no errors.  Check the log files
        # before the VCF files.  If there's an error, the VCF files
        # may not be created.
        # ##### ERROR -------------------------------------------------------
        # ##### ERROR A GATK RUNTIME ERROR has occurred (version 2.2-25-g2a68
        # ##### ERROR
        # ##### ERROR Please visit the wiki to see if this is a known problem
        # ##### ERROR If not, please post the error, with stack trace, to the
        # ##### ERROR Visit our website and forum for extensive documentation
        # ##### ERROR commonly asked questions http://www.broadinstitute.org/
        # ##### ERROR
        # ##### ERROR MESSAGE: java.lang.IllegalArgumentException: Comparison
        # ##### ERROR -------------------------------------------------------
        for i, x in enumerate(jobs):
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, raw_vcf_outfile, vcf_outfile, \
                log_outfile = x
            # Pull out the error lines.
            x = [x for x in open(log_outfile)]
            x = [x for x in x if x.startswith("##### ERROR")]
            x = "".join(x)
            msg = "MuTect error [%s]:\n%s\n%s" % (
                cancer_sample, commands[i], x)
            assert not x, msg

        # Make sure output VCF files exist.
        x = [x[6] for x in jobs]
        filelib.assert_exists_many(x)

        # Fix the files.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, raw_vcf_outfile, vcf_outfile, \
                log_outfile = x
            alignlib.clean_mutect_vcf(
                normal_bamfile, cancer_bamfile, normal_sample, cancer_sample,
                raw_vcf_outfile, vcf_outfile)
            
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "mutect.vcf"
