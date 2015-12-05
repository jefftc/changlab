from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents

        #in_filenames = filelib.list_files_in_path(
        #    bam_node.identifier, endswith=".bam", case_insensitive=True)
        in_filenames = module_utils.find_bam_files(bam_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # java -Xmx5g -jar /usr/local/bin/picard/picard.jar ReorderSam \
        #   I=<input.bam> O=<output.bam> REFERENCE=ucsc.hg19.fasta
        picard_jar = module_utils.find_picard_jar("picard")

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            x = in_filename, out_filename
            jobs.append(x)
        
        # Make a list of commands.
        sq = shell.quote
        commands = []
        for x in jobs:
            in_filename, out_filename = x

            x = [
                "java", "-Xmx5g",
                "-jar", sq(picard_jar),
                "ReorderSam",
                "I=%s" % sq(in_filename),
                "O=%s" % sq(out_filename),
                "REFERENCE=%s" % ref.fasta_file_full,
                ]
            x = " ".join(x)
            commands.append(x)
            
        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            in_filename, out_filename = x
            filelib.assert_exists_nz(out_filename)

    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFiles_sorted' + original_file
        #return filename
        return "bam"


