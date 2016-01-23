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

        vcf_node, ref_node = antecedents
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf")
        assert vcf_filenames, "No .vcf files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        
        jobs = []
        for in_filename in vcf_filenames:
            p, f = os.path.split(in_filename)
            f, exp = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.grp" % f)
            log_filename = os.path.join(out_path, "%s.log" % f)
            recal_filename = os.path.join(
                out_path, "%s.recalibrate_SNP.recal" % f)
            tranches_filename = os.path.join(
                out_path, "%s.recalibrate_SNP.tranches" % f)
            rscript_filename = os.path.join(
                out_path, "%s.recalibrate_SNP_plots.R" % f)
            assert in_filename != out_filename
            x = (in_filename, log_filename, recal_filename, tranches_filename,
                 rscript_filename)
            jobs.append(x)

        # -resource:dbsnp,known=true,training=false,truth=false,prior=6.0
        #    dbsnp_135.b37.vcf
        # -resource:hapmap,known=false,training=true,truth=true,prior=15.0
        #    hapmap_3.3.b37.sites.vcf
        # -resource:1000G,known=false,training=true,truth=false,prior=10.0
        #    1000G_phase1.snps.high_confidence.vcf
        # -resource:omni,known=false,training=true,truth=false,prior=12.0
        #    1000G_omni2.5.b37.sites.vcf
        known_sites = []
        x1 = module_utils.get_user_option(
            user_options, "vcf_recal_dbsnp", not_empty=True,
            check_file=True)
        x2 = module_utils.get_user_option(
            user_options, "vcf_recal_mills_indels", not_empty=True,
            check_file=True)
        x3 = module_utils.get_user_option(
            user_options, "vcf_recal_1kg_indels", not_empty=True,
            check_file=True)
        x4 = module_utils.get_user_option(
            user_options, "vcf_recal_omni", not_empty=True,
            check_file=True)
        y1 = "resource:dbsnp,known=true,training=false,truth=false,prior=6.0"
        y2 = "resource:hapmap,known=false,training=true,truth=true,prior=15.0"
        y3 = "resource:1000G,known=false,training=true,truth=false,prior=10.0"
        y4 = "resource:omni,known=false,training=true,truth=false,prior=12.0"
        known_sites = [(y1, x1), (y2, x2), (y3, x3), (y4, x4)]

        # Names of annotations to be used for annotations.
        AN = ["DP", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum",
              "InbreedingCoeff"]
        TRANCHE = ["100.0", "99.9", "99.0", "90.0"]
    
        # Make a list of commands.
        commands = []
        for x in jobs:
            (in_filename, log_filename, recal_filename, tranches_filename,
             rscript_filename) = x
            x1 = known_sites
            x2 = [("an", x) for x in AN]
            x3 = [("tranche", x) for x in TRANCHE]
            unhash = x1 + x2 + x3
            x = alignlib.make_GATK_command(
                T="VariantRecalibrator", R=ref.fasta_file_full,
                input=in_filename, mode="SNP", recalFile=recal_filename,
                tranchesFile=tranches_filename, rscriptFile=rscript_filename,
                _UNHASHABLE=unhash)
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "recal_vcf"

