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

        bam_path = in_data.identifier
        assert os.path.exists(bam_path)
        assert os.path.isdir(bam_path)
        filelib.safe_mkdir(out_path)

        # Find all the BAM files.
        bam_filenames = filelib.list_files_in_path(
            bam_path, endswith=".bam", case_insensitive=True)

        jobs = []  # list of in_filename, out_filename
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            assert not os.path.exists(out_filename)
            x = in_filename, out_filename
            jobs.append(x)

        # Symlink the BAM files to the output path.
        for x in jobs:
            in_filename, out_filename = x
            os.symlink(in_filename, out_filename)

        # Index each of the files.
        sq = parallel.quote
        samtools = filelib.which_assert(config.samtools)
        commands = []
        for x in jobs:
            in_filename, out_filename = x
            cmd = [
                sq(samtools),
                "index",
                sq(out_filename),
                ]
            x = " ".join(cmd)
            commands.append(x)
        
        parallel.pshell(commands, max_procs=num_cores, path=out_path)
    

    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFolder_' + original_file
        #return filename
        return "indexed.bam"

