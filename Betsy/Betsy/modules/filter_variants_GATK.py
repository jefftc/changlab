from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel

        vcf_node = in_data
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", not_empty=True)
        assert vcf_filenames, "No VCF files found."
        filelib.safe_mkdir(out_path)
        metadata = {}

        # Figure out whether the user wants SNPs or INDELs.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["snp", "indel"]
        metadata["filter"] = vartype

        jobs = []  # list of filelib.GenericObject
        for in_filename in vcf_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            x = filelib.GenericObject(
                in_filename=in_filename, out_filename=out_filename)
            jobs.append(x)

        # Filter each of the VCF files.
        jobs2 = []
        for j in jobs:
            args = vartype, j.in_filename, j.out_filename
            x = filter_by_vartype, args, {}
            jobs2.append(x)
        parallel.pyfun(jobs2, num_procs=num_cores)
        metadata["num_cores"] = num_cores

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "GATK.vcf"


def is_snp(var):
    # Doesn't seem to be annotated in the record.  Look at ref and alt
    # alleles.
    # REF  ALT
    # A    C     SNP
    # AG   A     INDEL
    assert len(var.ref) >= 1
    if len(var.ref) > 1:
        return False
    for x in var.alt:
        assert len(x) >= 1
        if len(x) > 1:
            return False
    return True


def is_indel(var):
    return not is_snp(var)


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
