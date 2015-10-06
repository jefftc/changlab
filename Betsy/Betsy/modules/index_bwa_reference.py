from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import config
        from genomicode import hashlib
        from Betsy import module_utils
        
        in_filename = in_data.identifier
        module_utils.safe_mkdir(out_path)
        bwa = module_utils.which_assert(config.bwa)

        # bwa index <out_stem.fa>
        # Makes files:
        # <out_stem>.fa.amb .ann .bwt .pac .sa

        out_stem = user_options.get("assembly", "genome")
        out_stem = hashlib.hash_var(out_stem)

        # Copy the in_filename to the out_path.
        assembly_filename = os.path.join(out_path, "%s.fa" % out_stem)
        shutil.copyfile(in_filename, assembly_filename)

        sq = module_utils.shellquote
        cwd = os.getcwd()
        try:
            os.chdir(out_path)

            cmd = [
                sq(bwa),
                "index",
                sq(assembly_filename),
                ]
            module_utils.run_single(cmd)
        finally:
            os.chdir(cwd)

        # Make sure the indexing worked properly.
        EXTENSIONS = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for ext in EXTENSIONS:
            f = "%s%s" % (assembly_filename, ext)
            assert module_utils.exists_nz(f), "Missing: %s" % f


    def name_outfile(self, antecedents, user_options):
        # Should name outfile based on the assembly.
        #from genomicode import hashlib
        #x = user_options.get("assembly", "genome")
        #x = hashlib.hash_var(x)
        #return "%s.bwa" % x
        return "reference.bwa"
