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

        bwa = filelib.which_assert(config.bwa)
        ref = alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)

        # bwa index <out_stem.fa>
        # Makes files:
        # <out_stem>.fa.amb .ann .bwt .pac .sa

        sq = parallel.quote
        cmd = [
            sq(bwa),
            "index",
            sq(ref.fasta_file_full),
            ]
        parallel.sshell(cmd, path=out_path)

        # Make sure the indexing worked properly.
        EXTENSIONS = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for ext in EXTENSIONS:
            f = "%s%s" % (ref.fasta_file_full, ext)
            assert filelib.exists_nz(f), "Missing: %s" % f

    def name_outfile(self, antecedents, user_options):
        return "reference.bwa"
