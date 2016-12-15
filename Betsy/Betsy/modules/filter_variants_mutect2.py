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
        return "mutect2.vcf"


def is_snp(var):
    from genomicode import vcflib
    # REF and ALT are one base each.
    # T  C
    # A  G
    assert type(var.ref) is type("")
    assert len(var.alt) == 1
    ref, alt = var.ref, var.alt[0]
    assert len(ref) >= 1 and len(alt) >= 1
    return len(ref) == 1 and len(alt) == 1


def is_indel(var):
    from genomicode import vcflib
    # REF and ALT are not one base
    # TG  T
    # CT  C
    # G   GCTATCAGTAAGTA
    assert type(var.ref) is type("")
    assert len(var.alt) == 1
    ref, alt = var.ref, var.alt[0]
    assert len(ref) >= 1 and len(alt) >= 1
    return len(ref) > 1 or len(alt) > 1


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
