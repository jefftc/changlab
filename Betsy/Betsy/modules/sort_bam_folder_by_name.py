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
        from genomicode import shell
        from Betsy import module_utils

        in_path = in_data.identifier
        filelib.safe_mkdir(out_path)

        in_filenames = filelib.list_files_in_path(
            in_path, endswith=".bam", case_insensitive=True)
        assert in_filenames, "No .bam files."

        # list of (in_filename, temp_prefix, log_filename, out_filename)
        jobs = []
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, f)
            temp_prefix = "temp_%s" % f
            log_filename = os.path.join(out_path, "%s.log" % s)
            assert in_filename != out_filename
            x = in_filename, temp_prefix, log_filename, out_filename
            jobs.append(x)

        # Make a list of samtools commands.
        samtools = filelib.which_assert(config.samtools)
        sq = shell.quote
        commands = []
        for x in jobs:
            in_filename, temp_prefix, log_filename, out_filename = x

            # samtools sort -n <in_filename> <out_filestem>
            # .bam automatically added to <out_filestem>, so don't
            # need it.
            x = out_filename
            assert x.endswith(".bam")
            x = x[:-4]
            out_filestem = x
            
            x = [
                samtools,
                "sort", "-n",
                "-O", "bam", 
                "-T", temp_prefix,
                "-o", sq(out_filename), 
                sq(in_filename),
                ]
            x = "%s >& %s" % (" ".join(x), sq(log_filename))
            commands.append(x)
            
        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        return "bam.sorted"

