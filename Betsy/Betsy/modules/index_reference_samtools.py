from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        from genomicode import config
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        
        samtools = filelib.which_assert(config.samtools)
        ref = alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)

        ## fa_filenames = module_utils.find_fasta_files(out_path)
        ## # Filter out the FASTA files created by RSEM indexing.
        ## # <assembly>.idx.fa
        ## # <assembly>.n2g.idx.fa
        ## # <assembly>.transcripts.fa
        ## # Could these end with ".fasta"?
        ## x = fa_filenames
        ## x = [x for x in x if not x.endswith(".idx.fa")]
        ## x = [x for x in x if not x.endswith(".n2g.idx.fa")]
        ## x = [x for x in x if not x.endswith(".transcripts.fa")]
        ## fa_filenames = x
        ## assert fa_filenames, "Could not find reference genome."
        ## assert len(fa_filenames) == 1, "Found multiple reference genomes."
        ## reference_filename = fa_filenames[0]

        # samtools faidx <ref>.fa
        # Makes files:
        # <ref>.fa.fai

        sq = parallel.quote
        cmd = [
            sq(samtools),
            "faidx",
            sq(ref.fasta_file_full),
            ]
        parallel.sshell(cmd, path=out_path)
        
        # Check to make sure index was created successfully.
        f = "%s.fai" % ref.fasta_file_full
        assert filelib.exists_nz(f)


    def name_outfile(self, antecedents, user_options):
        # Should name outfile based on the assembly.
        return "reference.samtools"
