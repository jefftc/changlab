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
        # TODO: Figure out version.

        # Figure out whether the user wants SNPs or INDELs.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["all", "snp", "indel"]

        # sample -> bam filename
        sample2bamfile = mlib.root2filename(bam_filenames)
        # Make sure files exist for all the samples.
        mlib.assert_normal_cancer_samples(nc_match, sample2bamfile)

        # list of (cancer_sample, normal_bamfile, tumor_bamfile, orig_outfile,
        #          fixed_outfile, filtered_outfile)
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_bamfile = sample2bamfile[normal_sample]
            cancer_bamfile = sample2bamfile[cancer_sample]
            path, sample, ext = mlib.splitpath(cancer_bamfile)
            orig_outfile = opj(out_path, "%s.raw" % sample)
            fix_outfile = opj(out_path, "%s.standard" % sample)
            filter_outfile = opj(out_path, "%s.vcf" % sample)
            x = cancer_sample, normal_bamfile, cancer_bamfile, \
                orig_outfile, fix_outfile, filter_outfile
            jobs.append(x)

        # python /usr/local/museq/classify.py \
        #   normal:test31/normal.bam tumour:test31/tumor.bam \
        #   reference:genomes/Broad.hg19/Homo_sapiens_assembly19.fa \
        #   model:/usr/local/museq/model_v4.1.2.npz \
        #   --config /usr/local/museq/metadata.config \
        #   -o test51.vcf
        opj = os.path.join
        museq = mlib.get_config("museq", assert_exists=True)
        classify_py = opj(museq, "classify.py")
        model_file = opj(museq, "model_v4.1.2.npz")
        config_file = opj(museq, "metadata.config")
        filelib.assert_exists_nz(classify_py)
        filelib.assert_exists_nz(model_file)
        filelib.assert_exists_nz(config_file)

        # museq's config file generates a broken VCF file.  Fix it.
        fixed_config_file = "fixed.config"
        fix_config_file(config_file, fixed_config_file)
  
        # Generate the commands.
        sq = mlib.sq
        commands = []
        for x in jobs:
            cancer_sample, normal_bamfile, cancer_bamfile, \
                           raw_outfile, fix_outfile, vcf_outfile = x

            x = [
                "python",    # should allow user to specify python
                sq(classify_py),
                sq("normal:%s" % normal_bamfile),
                sq("tumour:%s" % cancer_bamfile),
                sq("reference:%s" % ref.fasta_file_full),
                sq("model:%s" % model_file),
                "--config", sq(fixed_config_file),
                "-o", sq(raw_outfile),
                ]
            x = " ".join(map(str, x))
            commands.append(x)
        # Not sure how much RAM this takes.  On Thunderbolts test,
        # took < 1 Gb.
        nc = mlib.calc_max_procs_from_ram(5, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["num_cores"] = nc
        metadata["commands"] = commands

        # JointSNVMix produces non-standard VCF files.  Fix this so it
        # will work with other programs downstream.
        for x in jobs:
            cancer_sample, normal_bamfile, cancer_bamfile, \
                           raw_outfile, fix_outfile, vcf_outfile = x
            fix_vcf_file(cancer_sample, raw_outfile, fix_outfile)

        # Filter each of the VCF files.
        for x in jobs:
            cancer_sample, normal_bamfile, cancer_bamfile, \
                           raw_outfile, fix_outfile, vcf_outfile = x
            filter_by_vartype(vartype, fix_outfile, vcf_outfile)
        metadata["filter"] = vartype

        x = [x[-1] for x in jobs]
        filelib.assert_exists_many(x)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "jointsnvmix.vcf"


def is_snp(var):
    from genomicode import vcflib
    return vcflib.is_pass_filter(var, FILTER_doesnotcontain="INDL")


def is_indel(var):
    from genomicode import vcflib
    return vcflib.is_pass_filter(var, FILTER_doesnotcontain="PASS")


def filter_by_vartype(vartype, infile, outfile):
    # Filter out snps or indels.
    import shutil
    from genomicode import vcflib

    assert vartype in ["all", "snp", "indel"]

    if vartype == "all":
        shutil.copy2(infile, outfile)
        return
    vcf = vcflib.read(infile)
    fn = is_snp
    if vartype == "indel":
        fn = is_indel
    vcf = vcflib.select_variants(vcf, fn)
    vcflib.write(outfile, vcf)


def fix_vcf_file(sample, infile, outfile):
    # JointSNVMix produces VCF files that don't have FORMAT and
    # <SAMPLE> columns.  Add them.
    from genomicode import vcflib

    vcf = vcflib.read(infile)
    matrix = vcf.matrix

    genotype_names = ["DP", "RD", "AD", "FREQ"]

    # Get the calls for each variant.
    all_genotypes = []  # one for each variant
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        call = vcflib.get_call(var, None)
        geno_dict = {
            "DP" : call.total_reads,
            "RD" : call.num_ref,
            "AD" : call.num_alt,
            "FREQ" : call.vaf,
            }
        x = vcflib._format_genotype(genotype_names, geno_dict)
        all_genotypes.append(x)
    

    # Add FORMAT.
    FORMAT_STRING = ":".join(genotype_names)
    assert "FORMAT" not in matrix
    matrix.headers.append("FORMAT")
    matrix.headers_h.append("FORMAT")
    matrix.header2annots["FORMAT"] = [FORMAT_STRING] * matrix.num_annots()
    
    # Add the sample.
    assert not vcf.samples
    assert sample not in matrix
    matrix.headers.append(sample)
    matrix.headers_h.append(sample)
    matrix.header2annots[sample] = all_genotypes
    vcf.samples = [sample]

    # Add the proper header lines.
    lines = [
        '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Allelic depth for the ref allele in the tumor sample">',
        '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allelic depth for the alt allele in the tumor sample">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=FREQ,Number=1,Type=Integer,Description="Variant allele frequency">',
        ]
    matrix.headerlines.extend(lines)

    vcflib.write(outfile, vcf)


def fix_config_file(infile, outfile):
    # The config_file that JointSNVMix comes with generates a
    # broken VCF file.  bcftools merge gives error:
    # [W::vcf_parse] FILTER 'INDL' is not defined in the header
    #
    # To fix, need to add line:
    # ##FILTER=<ID=INDL,Description="indel">
    import shutil

    indl_line = '##FILTER=<ID=INDL,Description="indel">\n'
    lines = open(infile).readlines()
    
    # See if the line is already there.
    x = lines
    x = [x for x in x if x.startswith("##FILTER")]
    x = [x for x in x if x.find("INDL") >= 0]
    if x:
        # Line already there.  No need to fix.
        shutil.copy2(infile, outfile)
        return

    # Add to the end of the "##" lines.
    for i in range(len(lines)):
        if not lines[i].startswith("##"):
            break
    assert i < len(lines), "Missing non-## line"
    lines.insert(i, indl_line)
    open(outfile, 'w').writelines(lines)
