from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import alignlib
        from genomicode import parallel
        from Betsy import module_utils as mlib
        
        fastq_node, sample_node, strand_node, reference_node = antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        assert fastq_files, "I could not find any FASTQ files."
        ref = alignlib.create_reference_genome(reference_node.identifier)
        stranded = mlib.read_stranded(strand_node.identifier)
        filelib.safe_mkdir(out_path)
        
        metadata = {}
        metadata["tool"] = "RSEM %s" % alignlib.get_rsem_version()

        # Make a list of the jobs to run.
        jobs = []  # list of sample, pair1, pair2, log_filename
        for x in fastq_files:
            sample, pair1, pair2 = x
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, log_filename
            jobs.append(x)

        s2fprob = {
            "unstranded" : None,
            "firststrand" : 0.0,
            "secondstrand" : 1.0,
            }
        assert stranded.stranded in s2fprob, "Unknown stranded: %s" % \
               stranded.stranded
        forward_prob = s2fprob[stranded.stranded]

        sq = parallel.quote
        commands = []
        for x in jobs:
            sample, pair1, pair2, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_rsem_command(
                ref.fasta_file_full, sample, pair1, fastq_file2=pair2,
                forward_prob=forward_prob, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)
        metadata["commands"] = commands
        metadata["num cores"] = num_cores
        # Need to run in out_path.  Otherwise, files will be everywhere.
        parallel.pshell(commands, max_procs=num_cores, path=out_path)

        # Make sure the analysis completed successfully.
        x1 = [x[-1] for x in jobs]
        x2 = [os.path.join(out_path, "%s.genes.results" % x[0]) for x in jobs]
        x = x1 + x2
        filelib.assert_exists_nz_many(x)
        
        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "rsem"
