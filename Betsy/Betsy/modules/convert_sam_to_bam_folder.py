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
        from Betsy import module_utils

        filelib.safe_mkdir(out_path)

        in_path = module_utils.unzip_if_zip(in_data.identifier)
        x = filelib.list_files_in_path(in_path)
        x = [x for x in x if x.lower().endswith(".sam")]
        sam_filenames = x
        assert sam_filenames, "No .sam files."

        samtools = filelib.which_assert(config.samtools)

        jobs = []  # list of (sam_filename, bam_filename)
        for sam_filename in sam_filenames:
            p, f = os.path.split(sam_filename)
            assert f.endswith(".sam")
            f = f.replace(".sam", ".bam")
            bam_filename = os.path.join(out_path, f)
            x = sam_filename, bam_filename
            jobs.append(x)
        
        # Make a list of samtools commands.
        sq = parallel.quote
        commands = []
        for x in jobs:
            sam_filename, bam_filename = x

            # samtools view -bS -o <bam_filename> <sam_filename>
            x = [
                samtools,
                "view",
                "-bS",
                "-o", sq(bam_filename),
                sq(sam_filename),
                ]
            x = " ".join(x)
            commands.append(x)
            
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sam_filename, bam_filename = x
            assert filelib.exists_nz(bam_filename), \
                   "Missing: %s" % bam_filename
    

    def name_outfile(self, antecedents, user_options):
        return "bam"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #return 'bamFiles_' + original_file
