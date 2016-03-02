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
        from genomicode import hashlib
        from Betsy import module_utils
        
        bam_filenames = module_utils.find_bam_files(in_data.identifier)
        assert bam_filenames, "No .bam files."
        filelib.safe_mkdir(out_path)
        
        jobs = []  # list of (in_filename, sample, log_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            sample = hashlib.hash_var(s)
            log_filename = os.path.join(out_path, "%s.log" % s)
            out_filename = os.path.join(out_path, f)
            x = in_filename, sample, log_filename, out_filename
            jobs.append(x)
        
        gid = "group1"
        library = "library"
        platform_unit = "platform"
        #sample = "sample"
        platform = "illumina"

        # java -Xmx5g -jar AddOrReplaceReadGroups.jar
        #   I=<input.sam or .bam> O=<output.bam> ID=<group ID>
        #   LB=<group library> PU=<platform unit> SM=<group sample name>
        #   PL=<platform> CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
        picard_jar = alignlib.find_picard_jar("picard")

        # Make a list of commands.
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, sample, log_filename, out_filename = x
            x = [
                "java", "-Xmx5g",
                "-jar", sq(picard_jar),
                "AddOrReplaceReadGroups", 
                "I=%s" % sq(in_filename),
                "O=%s" % sq(out_filename),
                "ID=%s" % gid,
                "LB=%s" % library,
                "PU=%s" % platform_unit,
                "SM=%s" % sample,
                "PL=%s" % platform,
                #"CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        #for x in commands[:5]:
        #    print x
        #import sys; sys.exit(0)
            
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFiles_sorted' + original_file
        #return filename
        return "bam"


