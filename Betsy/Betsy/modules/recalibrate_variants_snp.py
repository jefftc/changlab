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

        vcf_node, ref_node, report_node = antecedents
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", not_empty=True)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        raise NotImplementedError

        recal_filenames = filelib.list_files_in_path(
            report_node.identifier, endswith="recalibrate_SNP.recal",
            not_empty=True)
        tranches_filenames = filelib.list_files_in_path(
            report_node.identifier, endswith="recalibrate_SNP.tranches",
            not_empty=True)
        assert len(vcf_filenames) == len(recal_filenames), \
               "Should have a recalibration file for each .vcf file."
        assert len(vcf_filenames) == len(tranches_filenames), \
               "Should have a tranches file for each .vcf file."
        
        sample2vcffilename = {}
        for filename in vcf_filenames:
            p, f = os.path.split(filename)
            sample, ext = os.path.splitext(f)
            assert sample not in sample2vcffilename
            sample2vcffilename[sample] = filename
        sample2recalfilename = {}
        for filename in recal_filenames:
            # <path>/<sample>.recalibrate_SNP.recal
            p, f = os.path.split(filename)
            n, ext = os.path.splitext(f)
            assert n.endswith(".recalibrate_SNP.recal")
            sample = n[:-len(".recalibrate_SNP.recal")]
            assert sample not in sample2recalfilename
            sample2recalfilename[sample] = filename
        sample2tranchesfilename = {}
        for filename in tranches_filenames:
            # <path>/<sample>.recalibrate_SNP.tranches
            p, f = os.path.split(filename)
            n, ext = os.path.splitext(f)
            assert n.endswith(".recalibrate_SNP.tranches")
            sample = n[:-len(".recalibrate_SNP.tranches")]
            assert sample not in sample2tranchesfilename
            sample2tranchesfilename[sample] = filename
        

        missing1 = [
            x for x in sample2vcffilename if x not in sample2recalfilename]
        missing2 = [
            x for x in sample2vcffilename if x not in sample2tranchesfilename]
        assert not missing1, \
               "Missing recalibration files for %d vcf files." % len(missing1)
        assert not missing2, \
               "Missing tranches files for %d vcf files." % len(missing2)

        # list of (in_filename, recal_filename, tranches_filename,
        #          log_filename, out_filename)
        jobs = []
        for sample in sample2vcffilename:
            in_filename = sample2vcffilename[sample]
            recal_filename = sample2recalfilename[sample]
            tranches_filename = sample2tranchesfilename[sample]
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, f)
            log_filename = os.path.join(out_path, "%s.log" % s)
            x = in_filename, recal_filename, tranches_filename, log_filename, \
                out_filename
            jobs.append(x)

        # java -jar GenomeAnalysisTK.jar \
        #   -T ApplyRecalibration \
        #   -mode SNP \
        #   --ts_filter_level 99.0 \
        #   -R reference.fa \
        #   -input raw_variants.vcf \
        #   -recalFile recalibrate_SNP.recal \
        #   -tranchesFile recalibrate_SNP.tranches \
        #   -o recalibrated_snps_raw_indels.vcf

        # Make a list of commands.
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, recal_filename, tranches_filename, log_filename, \
                         out_filename = x
            x = alignlib.make_GATK_command(
                T="ApplyRecalibration", mode="SNP",
                ts_filter_level="99.0", R=sq(ref.fasta_file_full),
                input=sq(in_filename), recalFile=sq(recal_filename),
                tranchesFile=sq(tranches_filename), o=sq(out_filename))
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "recalibrated.vcf"

