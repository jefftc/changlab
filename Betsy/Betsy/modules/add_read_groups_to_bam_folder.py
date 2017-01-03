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
        metadata = {}
        
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            sample = hashlib.hash_var(s)
            log_filename = os.path.join(out_path, "%s.log" % s)
            out_filename = os.path.join(out_path, f)
            x = filelib.GenericObject(
                in_filename=in_filename,
                sample=sample,
                log_filename=log_filename,
                out_filename=out_filename)
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
        for j in jobs:
            x = [
                "java", "-Xmx5g",
                "-jar", sq(picard_jar),
                "AddOrReplaceReadGroups", 
                "I=%s" % sq(j.in_filename),
                "O=%s" % sq(j.out_filename),
                "ID=%s" % gid,
                "LB=%s" % library,
                "PU=%s" % platform_unit,
                "SM=%s" % j.sample,
                "PL=%s" % platform,
                #"CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, sq(j.log_filename))
            commands.append(x)
            
        parallel.pshell(commands, max_procs=num_cores)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores

        # Make sure the analysis completed successfully.
        # Make sure outfiles exist.
        out_filenames = [j.out_filename for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

        # Check the log files to make sure there are no error.
        for j in jobs:
            check_log_file(j.log_filename)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "bam"


def check_log_file(filename):
    from genomicode import filelib
    # Log file format:
    # [Sat Dec 31 19:29:27 CST 2016] picard.sam.AddOrReplaceReadGroups INPUT=
    # [Sat Dec 31 19:29:27 CST 2016] Executing as jchang@velocitron.uth.tmc.e
    # INFO    2016-12-31 19:29:27     AddOrReplaceReadGroups  Created read gr
    # INFO    2016-12-31 19:29:42     AddOrReplaceReadGroups  Processed     1
    # INFO    2016-12-31 19:29:58     AddOrReplaceReadGroups  Processed     2
    # [...]
    # [Sat Dec 31 19:48:14 CST 2016] picard.sam.AddOrReplaceReadGroups done. 
    # Runtime.totalMemory()=1609564160
    #
    # Sometimes these lines are interspersed.  Probably OK not to flag.
    # Ignoring SAM validation error: ERROR: Read name HWI-ST1120:331:C6VW5ACX
    #
    # ERROR: Sometimes see exceptions.
    # Exception in thread "main" java.lang.RuntimeException: BGZF file has in
    #         at htsjdk.samtools.util.BlockGunzipper.unzipBlock(BlockGunzippe
    # [...]

    # The log file should not be empty.
    filelib.assert_exists_nz(filename)

    lines = open(filename).readlines()
    # Make sure there's no exception.
    i_exception = None
    for i in range(len(lines)):
        if lines[i].startswith("Exception in thread "):
            i_exception = i
            break
    if i_exception is None:
        return
    x = "".join(lines[i:]).strip()
    raise AssertionError, "Exception in Picard output:\n%s" % x
