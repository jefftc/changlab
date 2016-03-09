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

        in_path = in_data.identifier
        filelib.safe_mkdir(out_path)

        in_filenames = filelib.list_files_in_path(
            in_path, endswith=".bam", case_insensitive=True)
        assert in_filenames, "No .bam files."

        jobs = []  # list of (in_filename, tag_dir, log_file)
        seen = {}
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            x = f
            assert x.endswith(".bam")
            x = x[:-4]
            sample = x
            assert sample not in seen
            seen[sample] = 1

            tag_dir = sample
            log_file = "%s.log" % sample
            x = in_filename, tag_dir, log_file
            jobs.append(x)

        # Get the command.
        homer_path = filelib.which_assert(config.homer_path)
        x = os.path.join(homer_path, "bin", "makeTagDirectory")
        assert filelib.exists_nz(x)
        make_tag_directory = x
        
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, tag_dir, log_file = x

            # makeTagDirectory <tag_dir> <bam_file> >& log_file
            x = [
                sq(make_tag_directory),
                sq(tag_dir),
                sq(in_filename),
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, log_file)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores, path=out_path)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]
        x = [os.path.join(out_path, x) for x in x]
        log_filenames = x
        x = [x[1] for x in jobs]
        x = [os.path.join(out_path, x, "tagInfo.txt") for x in x]
        info_filenames = x
        x = log_filenames + info_filenames
        filelib.assert_exists_nz_many(x)

    
    def name_outfile(self, antecedents, user_options):
        return "tag_directories"
