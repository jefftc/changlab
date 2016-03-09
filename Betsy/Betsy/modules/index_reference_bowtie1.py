from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib

        ref = alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)
        bowtie_build = filelib.which_assert(config.bowtie_build)

        # bowtie-build <ref.fa> <name>
        # Makes files:
        # <name>.[1234].ebwt
        # <name>.rev.[12].ebwt

        sq = parallel.quote
        cmd = [
            sq(bowtie_build),
            sq(ref.fasta_file_full),
            ref.name,
            ]
        parallel.sshell(cmd, path=out_path)

        # Check to make sure index was created successfully.
        f = os.path.join(out_path, "%s.1.ebwt" % ref.name)
        assert filelib.exists_nz(f)


    def name_outfile(self, antecedents, user_options):
        return "reference.bowtie"
