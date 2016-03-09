from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # list of (in_filename, out_filename)
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, f)
            assert in_filename != out_filename
            x = in_filename, out_filename
            jobs.append(x)

        # Make a list of samtools commands.
        # Takes ~200 Mb per process, so should not be a big issue.
        samtools = filelib.which_assert(config.samtools)
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, out_filename = x

            # samtools calmd -b <in.bam> <ref.fasta> > <out.bam>
            x = [
                samtools,
                "calmd", "-b",
                sq(in_filename),
                sq(ref.fasta_file_full),
                ]
            x = "%s >& %s" % (" ".join(x), sq(out_filename))
            commands.append(x)
            
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)

    
    def name_outfile(self, antecedents, user_options):
        return "md.bam"


