from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        import call_variants_GATK

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

        jobs = []  # list of filelib.GenericObject
        for in_filename in vcf_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            x = filelib.GenericObject(
                in_filename=in_filename, out_filename=out_filename)
            jobs.append(x)

        # Filter each of the VCF files.
        for j in jobs:
            call_variants_GATK.filter_by_vartype(
                vartype, j.in_filename, j.out_filename)
        metadata["filter"] = vartype

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "platypus.vcf"
