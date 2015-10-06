from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from Betsy import module_utils

        in_file_or_path = in_data.identifier
        module_utils.copytree_or_file_into_tree(in_file_or_path, out_path)

        fa_filenames = module_utils.find_fasta_files(out_path)
        assert fa_filenames, "Could not find reference genome."
        assert len(fa_filenames) == 1, "Found multiple reference genomes."
        reference_filename = fa_filenames[0]

        f, e = os.path.splitext(reference_filename)
        out_filename = "%s.dict" % f

        # java -Xmx5g -jar /usr/local/bin/picard/CreateSequenceDictionary.jar \
        #   R=erdman.fa O=erdman.dict
        # <out>.dict

        picard_jar = module_utils.find_picard_jar("CreateSequenceDictionary")
        sq = module_utils.shellquote
        cmd = [
            "java", "-Xmx5g", "-jar", sq(picard_jar),
            # BUG: How to escape these filenames?
            "R=%s" % reference_filename,
            "O=%s" % out_filename,
            ]
        
        cwd = os.getcwd()
        try:
            os.chdir(out_path)
            module_utils.run_single(cmd)
        finally:
            os.chdir(cwd)

        # Make sure file was created successfully.
        assert module_utils.exists_nz(out_filename)


    def name_outfile(self, antecedents, user_options):
        # Should name outfile based on the assembly.
        return "reference.dict"
