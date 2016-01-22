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
        
        filelib.safe_mkdir(out_path)

        in_path = module_utils.unzip_if_zip(in_data.identifier)
        x = filelib.list_files_in_path(in_path)
        x = [x for x in x if x.lower().endswith(".bam")]
        in_filenames = x
        assert in_filenames, "No .bam files."

        # java -Xmx5g -jar MarkDuplicates.jar
        #   I=<input.sam or .bam> O=<output.bam> 
        #   METRICS_FILE=metricsFile CREATE_INDEX=true 
        #   VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
        #picard_jar = module_utils.find_picard_jar("MarkDuplicates")
        picard_jar = alignlib.find_picard_jar("picard")

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
                "MarkDuplicates",
                "I=%s" % sq(in_filename),
                "O=%s" % sq(out_filename),
                "METRICS_FILE=metricsFile",
                #"CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
                "REMOVE_DUPLICATES=true",
                ]
            x = " ".join(x)
            commands.append(x)
            
        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        return "bam"


