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
        from genomicode import parselib
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

        sample2bamfile = {}  # sample -> bam filename
        for filename in bam_filenames:
            path, sample, ext = mlib.splitpath(filename)
            sample2bamfile[sample] = filename

        # Make sure files exist for all the samples.
        all_samples = []
        for (normal_sample, cancer_sample) in nc_match:
            if normal_sample not in all_samples:
                all_samples.append(normal_sample)
            if cancer_sample not in all_samples:
                all_samples.append(cancer_sample)
        missing = [x for x in all_samples if x not in sample2bamfile]
        x = parselib.pretty_list(missing, max_items=5)
        assert not missing, "Missing BAM files for samples: %s" % x

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
            vcf_outfile = opj(out_path, "%s.vcf" % sample)
            log_outfile = opj(out_path, "%s.log" % sample)
            x = cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, vcf_outfile, log_outfile
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
            cancer_sample, normal_bamfile, cancer_bamfile, \
                call_outfile, cov_outfile, vcf_outfile, log_outfile = x

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
                vcf=sq(vcf_outfile),
                _UNHASHABLE=UNHASHABLE,
                )
            x = "%s >& %s" % (x, log_outfile)
            commands.append(x)
        nc = mlib.calc_max_procs_from_ram(15, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["num_cores"] = nc
        metadata["commands"] = commands
        
        x = [x[5] for x in jobs]
        filelib.assert_exists_many(x)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "mutect.vcf"
