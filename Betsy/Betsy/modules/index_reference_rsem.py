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
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        
        rsem_prepare = filelib.which_assert(config.rsem_prepare)
        ref = alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)

        gtf_file = user_options.get("gtf_file")
        assert gtf_file
        assert filelib.exists_nz(gtf_file), "File not found: %s" % gtf_file

        # rsem-prepare-reference --bowtie --bowtie2 --gtf gtf02.gtf
        #   <reference.fa> <reference_name>
        # <reference_name>.[1234].ebwt    # Bowtie1.
        # <reference_name>.rev.[12].ebwt
        # <reference_name>.[1234].bt2     # Bowtie2.
        # <reference_name>.rev.[12].bt2
        # <reference_name>.chrlist        # RSEM.
        # <reference_name>.grp
        # <reference_name>.idx.fa
        # <reference_name>.n2g.idx.fa
        # <reference_name>.seq
        # <reference_name>.ti
        # <reference_name>.transcripts.fa
        sq = parallel.quote
        # TODO: Need to test what happens if bowtie or bowtie2
        # can't be found.
        cmd = [
            sq(rsem_prepare),
            "--bowtie",
            "--bowtie2",
            "--gtf", sq(gtf_file),
            sq(ref.fasta_file_full),
            ref.name,
            ]
        parallel.sshell(cmd, path=out_path)

        # Copy the GTF file into the output path.
        shutil.copy2(gtf_file, out_path)

        assembly = ref.name
        # Check to make sure index was created successfully.
        x1 = ["%s.%d.ebwt" % (assembly, i+1) for i in range(4)]
        x2 = ["%s.rev.%d.ebwt" % (assembly, i+1) for i in range(2)]
        x3 = ["%s.%d.bt2" % (assembly, i+1) for i in range(4)]
        x4 = ["%s.rev.%d.bt2" % (assembly, i+1) for i in range(2)]
        x5 = [
            "%s.chrlist" % assembly,
            "%s.grp" % assembly,
            "%s.idx.fa" % assembly,
            "%s.n2g.idx.fa" % assembly,
            "%s.seq" % assembly,
            "%s.ti" % assembly,
            "%s.transcripts.fa" % assembly,
            ]
        x = x1 + x2 + x3 + x4 + x5
        index_files = [os.path.join(out_path, x) for x in x]
        filelib.assert_exists_nz_many(index_files)


    def name_outfile(self, antecedents, user_options):
        return "reference.rsem"
