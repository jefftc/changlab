from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        fastq_node, sample_node, reference_node = antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        ref = alignlib.create_reference_genome(reference_node.identifier)
        assert os.path.exists(ref.fasta_file_full)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "bowtie2 %s" % alignlib.get_bowtie2_version()

        # Make a list of the jobs to run.
        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            sam_filename = os.path.join(out_path, "%s.sam" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, sam_filename, log_filename
            jobs.append(x)
        
        sq = mlib.sq
        commands = []
        for x in jobs:
            sample, pair1, pair2, sam_filename, log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_bowtie2_command(
                ref.fasta_file_full, pair1, fastq_file2=pair2,
                sam_file=sam_filename, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-2] for x in jobs]
        filelib.assert_exists_nz_many(x)
            
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "alignments.bowtie2"
