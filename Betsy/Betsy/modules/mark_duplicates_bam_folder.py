from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib
        
        bam_filenames = mlib.find_bam_files(in_data.identifier)
        assert bam_filenames, "No .bam files."
        filelib.safe_mkdir(out_path)
        metadata = {}

        #in_path = module_utils.unzip_if_zip(in_data.identifier)
        #x = filelib.list_files_in_path(in_path)
        #x = [x for x in x if x.lower().endswith(".bam")]
        #in_filenames = x
        #assert in_filenames, "No .bam files."

        jobs = []  # list of (in_filename, log_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % s)
            out_filename = os.path.join(out_path, f)
            x = in_filename, log_filename, out_filename
            jobs.append(x)
            
        
        # java -Xmx5g -jar MarkDuplicates.jar
        #   I=<input.sam or .bam> O=<output.bam> 
        #   METRICS_FILE=metricsFile CREATE_INDEX=true 
        #   VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
        #picard_jar = module_utils.find_picard_jar("MarkDuplicates")
        picard_jar = alignlib.find_picard_jar("picard")

        # Make a list of commands.
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filename = x

            x = [
                "java", "-Xmx20g",
                "-jar", sq(picard_jar),
                "MarkDuplicates",
                "I=%s" % sq(in_filename),
                "O=%s" % sq(out_filename),
                "METRICS_FILE=metricsFile",
                #"CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
                #"REMOVE_DUPLICATES=true",
                "REMOVE_DUPLICATES=false",
                # BAM files should be sorted.
                # Actuallly, DEPRECATED now.
                #"ASSUME_SORTED=true",
                "TMP_DIR=%s" % sq(out_path),
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        # java may use additional threads for garbage collection.
        # https://sourceforge.net/p/picard/wiki/Main_Page/

        # Takes ~10 Gb per process (with -Xmx5g).  Increased the RAM
        # used to 20 Gb because 5 Gb was not enough for some files.
        nc = mlib.calc_max_procs_from_ram(25, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["commands"] = commands
        metadata["num_cores"] = nc

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "dups.bam"


