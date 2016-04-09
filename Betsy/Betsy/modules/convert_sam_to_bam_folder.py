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
        from Betsy import module_utils as mlib

        sam_filenames = mlib.find_sam_files(in_data.identifier)
        assert sam_filenames, "No .sam files."
        filelib.safe_mkdir(out_path)
        metadata = {}

        samtools = mlib.findbin("samtools")

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
                sq(samtools),
                "view",
                "-bS",
                "-o", sq(bam_filename),
                sq(sam_filename),
                ]
            x = " ".join(x)
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        return metadata
    

    def name_outfile(self, antecedents, user_options):
        return "bam"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #return 'bamFiles_' + original_file
