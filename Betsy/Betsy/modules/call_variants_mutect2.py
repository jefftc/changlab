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
        import call_somatic_varscan

        bam_node, nc_node, ref_node, interval_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.assert_exists_nz(interval_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        # TODO: Figure out GATK version.

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

        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            vcf_outfile = opj(out_path, "%s.vcf" % sample)
            log_outfile = opj(out_path, "%s.log" % sample)
            x = filelib.GenericObject(
                normal_sample=normal_sample,
                cancer_sample=cancer_sample,
                normal_bamfile=normal_bamfile,
                cancer_bamfile=cancer_bamfile,
                vcf_outfile=vcf_outfile,
                log_outfile=log_outfile)
            jobs.append(x)

        # java -jar GenomeAnalysisTK.jar \
        #   -T MuTect2 \
        #   -R reference.fasta \
        #   -I:tumor tumor.bam \
        #   -I:normal normal.bam \
        #   [--dbsnp dbSNP.vcf] \
        #   [--cosmic COSMIC.vcf] \
        #   [-L targets.interval_list] \
        #   -o output.vcf

        # Generate the commands.
        sq = mlib.sq
        commands = []
        for j in jobs:
            UNHASHABLE = [
                ("I:normal", sq(normal_bamfile)),
                ("I:tumor", sq(cancer_bamfile)),
                # --dbsnp and --cosmic use two dashes, for some
                # reason.  Since make_GATK_command only uses one dash,
                # add one manually.
                ("-dbsnp", sq(dbsnp_file)),
                ("-cosmic", sq(cosmic_file)),
                ]
            x = alignlib.make_GATK_command(
                T="MuTect2",
                R=sq(ref.fasta_file_full),
                L=sq(interval_node.identifier),
                o=sq(j.vcf_outfile),
                _UNHASHABLE=UNHASHABLE,
                )
            x = "%s >& %s" % (x, j.log_outfile)
            commands.append(x)
        assert len(commands) == len(jobs)
        
        nc = mlib.calc_max_procs_from_ram(25, upper_max=num_cores)
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
        for i, j in enumerate(jobs):
            # Pull out the error lines.
            x = [x for x in open(j.log_outfile)]
            x = [x for x in x if x.startswith("##### ERROR")]
            x = "".join(x)
            msg = "MuTect2 error [%s]:\n%s\n%s" % (
                cancer_sample, commands[i], x)
            assert not x, msg

        # Make sure output VCF files exist.
        x = [x.vcf_outfile for x in jobs]
        filelib.assert_exists_many(x)

        # Mutect2 names the samples "NORMAL" and "TUMOR".  Replace
        # them with the actual names.
        for j in jobs:
            call_somatic_varscan._fix_normal_cancer_names(
                j.vcf_outfile, j.normal_sample, j.cancer_sample)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "mutect2.vcf"
