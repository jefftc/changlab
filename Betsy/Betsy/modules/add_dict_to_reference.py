from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        #picard_jar = module_utils.find_picard_jar("CreateSequenceDictionary")
        picard_jar = module_utils.find_picard_jar("picard")
        ref = alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)
        out_file = "%s.dict" % ref.fasta_file

        # java -Xmx5g -jar /usr/local/bin/picard/CreateSequenceDictionary.jar \
        #   R=erdman.fa O=erdman.dict
        # <out>.dict

        sq = shell.quote
        cmd = [
            "java", "-Xmx5g", "-jar", sq(picard_jar),
            "CreateSequenceDictionary",
            # BUG: How to escape these filenames?
            "R=%s" % ref.fasta_file_full,
            "O=%s" % out_file,
            ]
        shell.single(cmd, path=out_path)

        # Make sure file was created successfully.
        x = os.path.join(out_path, out_file)
        assert filelib.exists_nz(x), "File not found: %s" % out_file


    def name_outfile(self, antecedents, user_options):
        return "ref.dict"
