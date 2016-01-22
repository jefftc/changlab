from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (bam_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            f, exp = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.grp" % f)
            assert in_filename != out_filename
            x = in_filename, out_filename
            jobs.append(x)

        known_sites = []
        x1 = module_utils.get_user_option(
            user_options, "recal_known_sites1", not_empty=True,
            check_file=True)
        x2 = module_utils.get_user_option(
            user_options, "recal_known_sites2", check_file=True)
        x3 = module_utils.get_user_option(
            user_options, "recal_known_sites3", check_file=True)
        x = [x1, x2, x3]
        x = [x for x in x if x]
        known_sites = x
        assert known_sites
    
        # java -Xmx4g -jar GenomeAnalysisTK.jar
        #   -T BaseRecalibrator -R /Path/hg19.fa
        #    -knownSites /Path/bundle-1.5/hg19/dbsnp_135.hg19.vcf
        #    -I /Path/mySample.bam -o /Path/mySample_CovarTable_Recal.grp
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, out_filename = x
            x = [("knownSites", x) for x in known_sites]
            x = alignlib.make_GATK_command(
                T="BaseRecalibrator", R=ref.fasta_file_full,
                I=in_filename, o=out_filename, _UNHASHABLE=x)
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "recalibration.grp"

