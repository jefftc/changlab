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
        from Betsy import module_utils
        
        in_filenames = module_utils.find_bam_files(in_data.identifier)
        assert in_filenames, "No .bam files."
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "samtools %s" % alignlib.get_samtools_version()
        
        jobs = []  # list of (in_filename, temp_prefix, out_filename)
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            temp_prefix = "temp_%s" % f
            x = in_filename, temp_prefix, out_filename
            jobs.append(x)
        
        samtools = filelib.which_assert(config.samtools)

        # Calculate the number of threads per process.
        nc = module_utils.calc_max_procs_from_ram(4, upper_max=num_cores)
        num_threads = max(nc/len(jobs), 1)

        # Make a list of samtools commands.
        # Without -m, takes ~1 Gb per process.
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, temp_prefix, out_filename = x

            # Usage has changed.  Below no longer valid.
            # samtools sort <in_filename> <out_filestem>
            # .bam automatically added to <out_filestem>, so don't
            # need it.
            #x = out_filename
            #assert x.endswith(".bam")
            #x = x[:-4]
            #out_filestem = x
            
            x = [
                sq(samtools),
                "sort",
                "-O", "bam",
                "-T", sq(temp_prefix),
                "-m", "4G",    # Crashing, so try increasing memory.
                sq(in_filename),
                "-o", sq(out_filename),
                ]
            if num_threads > 1:
                x += ["-@", num_threads]
            x = " ".join(map(str, x))
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = nc

        parallel.pshell(commands, max_procs=nc)
        #for cmd in commands:
        #    parallel.sshell(cmd)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)
        
        return metadata

    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFiles_sorted' + original_file
        #return filename
        return "sorted.coord.bam"


