from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import parallel
        from genomicode import alignlib
        from genomicode import filelib
        from Betsy import module_utils as mlib
        
        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        metadata["tool"] = "samtools %s" % alignlib.get_samtools_version()

        # list of (in_filename, err_filename, out_filename)
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            sample, ext = os.path.splitext(f)
            err_filename = os.path.join(out_path, "%s.log" % sample)
            out_filename = os.path.join(out_path, "%s.pileup" % sample)
            x = in_filename, err_filename, out_filename
            jobs.append(x)

        # samtools mpileup -f [reference sequence] [BAM file(s)]
        #   > myData.mpileup
        samtools = mlib.findbin("samtools")
        sq = mlib.sq
        commands = []
        for x in jobs:
            in_filename, err_filename, out_filename = x
            
            x = [
                sq(samtools),
                "mpileup",
                "-f", sq(ref.fasta_file_full),
                ]
            x.append(sq(in_filename))
            x = " ".join(map(str, x))
            x = "%s 2> %s 1> %s" % (x, err_filename, out_filename)
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores)
        metadata["num_cores"] = num_cores
        metadata["commands"] = commands

        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)

        return metadata
        

    def name_outfile(self, antecedents, user_options):
        return "summary.pileup"
