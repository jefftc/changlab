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
        #from genomicode import hashlib
        from Betsy import module_utils
        
        in_filenames = module_utils.find_bam_files(in_data.identifier)
        assert in_filenames, "No .bam files."
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "samtools %s" % alignlib.get_samtools_version()

        jobs = []
        #seen = {}
        for i, in_filename in enumerate(in_filenames):
            p, f = os.path.split(in_filename)
            temp_prefix = "temp_%s" % f
            #temp_prefix = "temp_%s" % hashlib.hash_var(f)
            # Make sure no duplicates.
            #assert temp_prefix not in seen
            #seen[temp_prefix] = 1
            #temp_outfilename = "%d.bam" % i
            out_filename = os.path.join(out_path, f)
            x = filelib.GenericObject(
                in_filename=in_filename,
                temp_prefix=temp_prefix,
                #temp_outfilename=temp_outfilename,
                out_filename=out_filename)
            jobs.append(x)
            
        samtools = filelib.which_assert(config.samtools)

        # Calculate the number of threads per process.
        nc = module_utils.calc_max_procs_from_ram(4, upper_max=num_cores)
        num_threads = max(nc/len(jobs), 1)

        # Make a list of samtools commands.
        # Without -m, takes ~1 Gb per process.
        sq = parallel.quote
        commands = []
        for j in jobs:
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
                "-T", sq(j.temp_prefix),
                "-m", "4G",    # Crashing, so try increasing memory.
                sq(j.in_filename),
                #"-o", sq(j.temp_outfilename),
                "-o", sq(j.out_filename),
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

        #for j in jobs:
        #    # Move the temporary files to the final location.
        #    shutil.move(j.temp_outfilename, j.out_filename)

        # Make sure the analysis completed successfully.
        x = [j.out_filename for j in jobs]
        filelib.assert_exists_nz_many(x)
        
        return metadata

    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFiles_sorted' + original_file
        #return filename
        return "sorted.coord.bam"


