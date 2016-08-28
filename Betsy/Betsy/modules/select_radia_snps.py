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
        vcf_files = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", case_insensitive=True)
        filelib.safe_mkdir(out_path)
        metadata = {}

        jobs = []  # in_vcf_filename, out_vcf_filename
        for vcf_file in vcf_files:
            path, file_ = os.path.split(vcf_file)
            out_vcf_file = os.path.join(out_path, file_)
            x = vcf_file, out_vcf_file
            jobs.append(x)

        # Figure out whether the user wants SNPs or INDELs.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["all", "snp", "indel"]

        # Generate the commands.
        commands = []
        for x in jobs:
            in_vcf_file, out_vcf_file = x

            args = vartype, in_vcf_file, out_vcf_file
            x = filter_by_vartype, args, {}
            commands.append(x)
        parallel.pyfun(commands, num_procs=num_cores)
        metadata["num_cores"] = num_cores

        x = [x[-1] for x in jobs]
        filelib.assert_exists_many(x)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "radia.vcf"


def is_snp(var):
    # Want the records with INFO:VT == "SNP".
    vt = var.infodict.get("VT")
    assert vt in ["SNP", "INS", "DEL"]
    return vt == "SNP"


def is_indel(var):
    # Want the records with INFO:VT == "INS" or "DEL"
    vt = var.infodict.get("VT")
    assert vt in ["SNP", "INS", "DEL"]
    return vt in ["INS", "DEL"]


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
