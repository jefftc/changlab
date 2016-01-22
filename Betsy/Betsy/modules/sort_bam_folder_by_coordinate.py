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
        
        filelib.safe_mkdir(out_path)

        in_path = module_utils.unzip_if_zip(in_data.identifier)
        x = filelib.list_files_in_path(in_path)
        x = [x for x in x if x.lower().endswith(".bam")]
        in_filenames = x
        assert in_filenames, "No .bam files."

        samtools = filelib.which_assert(config.samtools)

        jobs = []  # list of (in_filename, temp_prefix, out_filename)
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            temp_prefix = "temp_%s" % f
            x = in_filename, temp_prefix, out_filename
            jobs.append(x)
        
        # Make a list of samtools commands.
        sq = shell.quote
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
                samtools,
                "sort",
                "-O", "bam",
                "-T", temp_prefix,
                sq(in_filename),
                "-o", sq(out_filename),
                ]
            x = " ".join(x)
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFiles_sorted' + original_file
        #return filename
        return "bam"


