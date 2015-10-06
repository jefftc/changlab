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
        
        in_file_or_path = in_data.identifier
        module_utils.copytree_or_file_into_tree(in_file_or_path, out_path)

        fa_filenames = module_utils.find_fasta_files(out_path)
        # Filter out the FASTA files created by RSEM indexing.
        # <assembly>.idx.fa
        # <assembly>.n2g.idx.fa
        # <assembly>.transcripts.fa
        # Could these end with ".fasta"?
        x = fa_filenames
        x = [x for x in x if not x.endswith(".idx.fa")]
        x = [x for x in x if not x.endswith(".n2g.idx.fa")]
        x = [x for x in x if not x.endswith(".transcripts.fa")]
        fa_filenames = x
        assert fa_filenames, "Could not find reference genome."
        assert len(fa_filenames) == 1, "Found multiple reference genomes."
        reference_filename = fa_filenames[0]

        # samtools faidx <ref>.fa
        # Makes files:
        # <ref>.fa.fai

        sq = module_utils.shellquote
        samtools = module_utils.which_assert(config.samtools)
        cmd = [
            sq(samtools),
            "faidx",
            sq(reference_filename),
            ]
        
        cwd = os.getcwd()
        try:
            os.chdir(out_path)
            module_utils.run_single(cmd)
        finally:
            os.chdir(cwd)

        # Check to make sure index was created successfully.
        f = "%s.fai" % reference_filename
        assert module_utils.exists_nz(f)


    def name_outfile(self, antecedents, user_options):
        # Should name outfile based on the assembly.
        return "reference.samtools"
