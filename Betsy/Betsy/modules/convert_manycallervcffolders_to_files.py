from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from Betsy import module_utils as mlib
        import merge_vcf_folder

        vcffolders_node = antecedents
        filelib.safe_mkdir(out_path)
        metadata = {}

        x = os.listdir(vcffolders_node.identifier)
        x = [x for x in x if x.endswith(".vcf")]
        assert x, "No VCF folders found: %s" % vcffolders_node.identifier
        x = [os.path.join(vcffolders_node.identifier, x) for x in x]
        vcf_folders = x

        jobs = []
        for folder in vcf_folders:
            path, root, ext = mlib.splitpath(folder)
            assert ext == ".vcf"
            caller = root
            vcf_filenames = filelib.list_files_in_path(
                folder, endswith=".vcf", toplevel_only=True)
            assert vcf_filenames, "No .vcf files: %s" % folder
            out_filename = os.path.join(out_path, "%s.vcf" % root)
            tmp_path = "%s.indexed.vcf" % caller
            x = filelib.GenericObject(
                caller=caller,
                vcf_filenames=vcf_filenames,
                out_filename=out_filename,
                tmp_path=tmp_path
                )
            jobs.append(x)

        for j in jobs:
            m = merge_vcf_folder.merge_vcf_files(
                j.vcf_filenames, j.out_filename, num_cores, j.tmp_path)
            if "commands" not in metadata:
                metadata["commands"] = []
            metadata["commands"].extend(m["commands"])

        x = [x.out_filename for x in jobs]
        filelib.assert_exists_many(x)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "multisample.vcf"
