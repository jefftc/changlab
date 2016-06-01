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

        bam_node, nc_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        # TODO: Figure out version.

        # sample -> bam filename
        sample2bamfile = mlib.root2filename(bam_filenames)
        # Make sure files exist for all the samples.
        mlib.assert_normal_cancer_samples(nc_match, sample2bamfile)

        # list of (normal_sample, cancer_sample, normal_bamfile, tumor_bamfile,
        #          vcf_outfile)
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            vcf_outfile = opj(out_path, "%s.vcf" % sample)
            x = normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                vcf_outfile
            jobs.append(x)

        # bam-somaticsniper -q 1 -Q 15 -G -L -F vcf \
        #   -f genomes/Broad.hg19/Homo_sapiens_assembly19.fa \
        #   test31/tumor.bam test31/normal.bam test41.vcf
        somaticsniper = mlib.get_config(
            "somaticsniper", which_assert_file=True)
  
        # Generate the commands.
        sq = mlib.sq
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           vcf_outfile = x

            x = [
                sq(somaticsniper),
                "-q", 1,
                "-Q", 15,
                "-G",
                "-L",
                "-F", "vcf",
                "-f", sq(ref.fasta_file_full),
                sq(cancer_bamfile),
                sq(normal_bamfile),
                sq(vcf_outfile),
                ]
            x = " ".join(map(str, x))
            commands.append(x)
        # Not sure how much RAM this takes.
        nc = mlib.calc_max_procs_from_ram(15, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["num_cores"] = nc
        metadata["commands"] = commands

        # SomaticSniper names the samples "NORMAL" and "TUMOR".
        # Replace them with the actual names.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           vcf_outfile = x
            call_somatic_varscan._fix_normal_cancer_names(
                vcf_outfile, normal_sample, cancer_sample)
        
        x = [x[-1] for x in jobs]
        filelib.assert_exists_many(x)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "somaticsniper.vcf"
