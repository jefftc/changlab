from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from Betsy import module_utils
        
        in_filename = in_data.identifier
        module_utils.safe_mkdir(out_path)
        bowtie2_build = module_utils.which_assert(config.bowtie2_build)
        

        # bowtie2-build <ref.fa> <output_stem>
        # Makes files:
        # <output_stem>.[1234].bt2
        # <output_stem>.rev.[12].bt2

        out_stem = user_options.get("assembly", "genome")

        # Figure out the output stem.
        # Not good, because this is often a weird betsy-defined name.
        #in_path, in_file = os.path.split(in_filename)
        #x = in_file
        #if x.lower().endswith(".fa"):
        #    x = x[:-3]
        #if x.lower().endswith(".fasta"):
        #    x = x[:-6]
        #out_stem = x

        cwd = os.getcwd()
        try:
            os.chdir(out_path)

            cmd = [
                bowtie2_build,
                in_filename,
                out_stem,
                ]
            module_utils.run_single(cmd)
        finally:
            os.chdir(cwd)

        # Check to make sure index was created successfully.
        f = os.path.join(out_path, "%s.1.bt2" % out_stem)
        assert module_utils.exists_nz(f)


    def name_outfile(self, antecedents, user_options):
        return "reference.bowtie2"
