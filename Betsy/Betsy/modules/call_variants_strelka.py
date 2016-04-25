from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
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
        # TODO: Figure out Strelka version.

        skip_depth_filter = False
        x = mlib.get_user_option(
            user_options, "strelka_skip_depth_filter",
            allowed_values=["no", "yes"], not_empty=True)
        if x == "yes":
            skip_depth_filter = True
        assert "vartype" in out_attributes, "Missing attribute: vartype"
        x = out_attributes["vartype"]
        assert x in ["snp", "indel"]
        vartype = x

        # sample -> bam filename
        sample2bamfile = mlib.root2filename(bam_filenames)
        # Make sure files exist for all the samples.
        mlib.assert_normal_cancer_samples(nc_match, sample2bamfile)

        # list of (normal_sample, cancer_sample, normal_bamfile, tumor_bamfile,
        #          config_file, output_dir
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            config_file = opj(out_path, "config.%s.ini" % cancer_sample)
            analysis_path = opj(out_path, "analysis.%s" % cancer_sample)
            x = normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                config_file, analysis_path
            jobs.append(x)

        # Make each of the config files.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           config_file, analysis_path = x
            _make_config_file(config_file, skip_depth_filter=skip_depth_filter)

        # Make the analysis directories.
        jobs2 = []
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           config_file, analysis_path = x
            fn = _make_analysis_directory
            args = (
                analysis_path, config_file, ref.fasta_file_full,
                normal_bamfile, cancer_bamfile)
            keywds = None
            jobs2.append((fn, args, keywds))
        parallel.pyfun(jobs2, num_procs=num_cores)

        # Run the analysis.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           config_file, analysis_path = x
            cmd = "make -j %d" % num_cores
            parallel.sshell(cmd, path=analysis_path)
        metadata["num_cores"] = num_cores

        # Make sure files exists.
        x = [x[-1] for x in jobs]
        x = [os.path.join(x, "results", "all.somatic.snvs.vcf") for x in x]
        filelib.assert_exists_nz_many(x)

        # Clean the VCF files and save into the out_path.
        for x in jobs:
            normal_sample, cancer_sample, normal_bamfile, cancer_bamfile, \
                           config_file, analysis_path = x
            # <analysis_path>/results/all.somatic.snvs.vcf
            # <analysis_path>/results/all.somatic.indels.vcf
            vartype2file = {
                "snp" : "all.somatic.snvs.vcf",
                "indel" : "all.somatic.indels.vcf",
                }
            assert vartype in vartype2file
            x = vartype2file[vartype]
            src_file = os.path.join(analysis_path, "results", x)
            dst_file = os.path.join(out_path, "%s.vcf" % cancer_sample)
            alignlib.clean_strelka_vcf(
                normal_sample, cancer_sample, src_file, dst_file)


        #metadata["commands"] = commands
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "strelka.vcf"


def _make_config_file(config_filename, skip_depth_filter=False):
    import os
    from genomicode import filelib
    from Betsy import module_utils as mlib

    strelka_path = mlib.get_config("strelka", assert_exists=True)

    src_config = os.path.join(
        strelka_path, "etc", "strelka_config_bwa_default.ini")
    filelib.exists_nz(src_config)
    lines = open(src_config).readlines()
    assert lines

    # Edit configure options.
    for i in range(len(lines)):
        x = lines[i]
        x = x.strip()
        line = x
    
        # Make sure skip_depth_filter is correct.
        # isSkipDepthFilters should be set to 1 to skip depth
        # filtration for whole exome or other targeted sequencing data
        # 
        # sSkipDepthFilters = 0
        if line.startswith("isSkipDepthFilters"):
            # isSkipDepthFilters = 0
            x = line.split()
            assert len(x) == 3
            assert x[1] == "="
            assert x[2] in ["0", "1"]
            if skip_depth_filter:
                x[2] = "1"
            else:
                x[2] = "0"
            line = " ".join(x)
        lines[i] = line

    lines = [x+"\n" for x in lines]  # replace newline that was stripped.
    open(config_filename, 'w').writelines(lines)


def _make_analysis_directory(analysis_path,
                             config_file, reference_fa, normal_bam, tumor_bam):
    import os
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_utils as mlib

    filelib.assert_exists_nz(config_file)
    filelib.assert_exists_nz(reference_fa)
    filelib.assert_exists_nz(normal_bam)
    filelib.assert_exists_nz(tumor_bam)

    strelka_path = mlib.get_config("strelka", assert_exists=True)
    config_pl = os.path.join(
        strelka_path, "bin", "configureStrelkaWorkflow.pl")
    filelib.assert_exists_nz(config_pl)

    # $STRELKA/bin/configureStrelkaWorkflow.pl \
    #   --normal=../test31.bam --tumor=../test32.bam \
    #   --ref=../genomes/Broad.hg19/Homo_sapiens_assembly19.fa \
    #   --config=./config.ini --output-dir=./myAnalysis
    sq = mlib.sq
    cmd = [
        sq(config_pl),
        "--normal", sq(normal_bam),
        "--tumor", sq(tumor_bam),
        "--ref", sq(reference_fa),
        "--config", sq(config_file),
        "--output-dir", sq(analysis_path),
        ]
    cmd = " ".join(cmd)
    parallel.sshell(cmd)


