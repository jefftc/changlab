from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # list of (in_filename, log_filename, err_filename, out_filename)
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            sample, ext = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            err_filename = os.path.join(out_path, "%s.err" % sample)
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = in_filename, log_filename, err_filename, out_filename
            jobs.append(x)

        buffer_size = 100000
        max_reads = 5E6
        # Running into errors sometimes, so increase these numbers.
        #   WARNING - Too many reads (5000000) in region
        #   1:500000-600000. Quitting now. Either reduce --bufferSize or
        #   increase --maxReads.
        buffer_size = buffer_size * 10
        max_reads = max_reads * 10
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, log_filename, err_filename, out_filename = x
            #nc = max(1, num_cores/len(jobs))
            x = alignlib.make_platypus_command(
                bam_file=in_filename,
                ref_file=ref.fasta_file_full,
                log_file=log_filename,
                out_file=out_filename,
                buffer_size=buffer_size, max_reads=max_reads)
            x = "%s >& %s" % (x, err_filename)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.  If not, try
        # to diagnose.
        for x in jobs:
            x, x, err_filename, out_filename = x
            if filelib.exists_nz(out_filename):
                continue
            for line in open(err_filename):
                if line.find("WARNING - Too many reads") >= 0:
                    print line,
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "platypus.vcf"
