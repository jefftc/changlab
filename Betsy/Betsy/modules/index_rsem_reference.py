from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import hashlib
        from Betsy import module_utils
        
        in_filename = in_data.identifier
        assert os.path.exists(in_filename)
        module_utils.safe_mkdir(out_path)
        rsem_prepare = module_utils.which_assert(config.rsem_prepare)

        gtf_file = user_options.get("gtf_file")
        assert gtf_file
        assert module_utils.exists_nz(gtf_file)

        assembly = user_options["assembly"]
        assert assembly
        # Make sure this can be a valid filename.
        assembly = hashlib.hash_var(assembly)

        # rsem-prepare-reference --bowtie --bowtie2 --gtf gtf02.gtf
        #   <reference.fa> <reference_name>
        # <reference_name>.[1234].ebwt    # Bowtie1
        # <reference_name>.rev.[12].ebwt
        # <reference_name>.[1234].bt2      # Bowtie2
        # <reference_name>.rev.[12].bt2
        # <reference_name>.chrlist
        # <reference_name>.grp
        # <reference_name>.idx.fa
        # <reference_name>.n2g.idx.fa
        # <reference_name>.seq
        # <reference_name>.ti
        # <reference_name>.transcripts.fa

        sq = module_utils.shellquote
        cwd = os.getcwd()
        try:
            os.chdir(out_path)

            # TODO: Need to test what happens if bowtie or bowtie2
            # can't be found.
            cmd = [
                sq(rsem_prepare),
                "--bowtie",
                "--bowtie2",
                "--gtf", sq(gtf_file),
                sq(in_filename),
                assembly,
                ]
            cmd = " ".join(cmd)
            module_utils.run_single(cmd)
        finally:
            os.chdir(cwd)

        # Check to make sure index was created successfully.
        index_files = []
        x1 = ["%s.%d.ebwt" % (assembly, i+1) for i in range(4)]
        x2 = ["%s.rev.%d.ebwt" % (assembly, i+1) for i in range(2)]
        x3 = ["%s.%d.bt2" % (assembly, i+1) for i in range(4)]
        x4 = ["%s.rev.%d.btw" % (assembly, i+1) for i in range(2)]
        index_files += x1 + x2 + x3 + x4
        index_files.append("%s.chrlist" % assembly)
        index_files.append("%s.grp" % assembly)
        index_files.append("%s.idx.fa" % assembly)
        index_files.append("%s.n2g.idx.fa" % assembly)
        index_files.append("%s.seq" % assembly)
        index_files.append("%s.ti" % assembly)
        index_files.append("%s.transcripts.fa" % assembly)
        for f in index_files:
            filename = os.path.join(out_path, f)
            assert module_utils.exists_nz(filename)


    def name_outfile(self, antecedents, user_options):
        return "reference.rsem"
