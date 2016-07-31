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

        # Figure out whether to align to genome or transcriptome.
        x = out_attributes["align_to"]
        assert x in ["genome", "transcriptome"]
        align_to_genome = (x == "genome")

        # Make a list of the jobs to run.
        jobs = []  # list of sample, pair1, pair2, log_filename
        for x in fastq_files:
            sample, pair1, pair2 = x
            results_filename = os.path.join(
                out_path, "%s.genes.results" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = filelib.GenericObject(
                sample=sample, pair1=pair1, pair2=pair2,
                results_filename=results_filename, 
                log_filename=log_filename)
            jobs.append(x)

        s2fprob = {
            "unstranded" : None,
            "firststrand" : 0.0,
            "secondstrand" : 1.0,
            }
        assert stranded.stranded in s2fprob, "Unknown stranded: %s" % \
               stranded.stranded
        forward_prob = s2fprob[stranded.stranded]

        # How much memory for bowtie.  May need to increase this if
        # there are lots of memory warnings in the log files:
        #   Warning: Exhausted best-first chunk memory for read
        #   ST-J00106:110:H5NY5BBXX:6:1101:18203:44675 1:N:0:1/1
        #   (patid 2076693); skipping read
        # Default is 64.
        # Seems like too high a value can cause problems.
        #chunkmbs = 4*1024   # Generates warnings.
        chunkmbs = 512

        # Get lots of warnings with bowtie:
        # Warning: Detected a read pair whose two mates have different names

        # Use STAR aligner instead.

        use_STAR = True

        sq = parallel.quote
        commands = []
        for j in jobs:
            # Debug: If the results file exists, don't run it again.
            if filelib.exists_nz(j.results_filename):
                continue
            # If using the STAR aligner, then most memory efficient
            # way is to let STAR take care of the multiprocessing.
            nc = max(1, num_cores/len(jobs))
            if use_STAR:
                nc = num_cores

            keywds = {}
            if use_STAR:
                keywds["align_with_star"] = True
            else:
                keywds["align_with_bowtie2"] = True
            x = alignlib.make_rsem_command(
                ref.fasta_file_full, j.sample, j.pair1, fastq_file2=j.pair2,
                forward_prob=forward_prob, output_genome_bam=align_to_genome,
                bowtie_chunkmbs=chunkmbs, num_threads=nc, **keywds)
            x = "%s >& %s" % (x, sq(j.log_filename))
            commands.append(x)
        metadata["commands"] = commands
        metadata["num cores"] = num_cores
        # Need to run in out_path.  Otherwise, files will be everywhere.
        nc = num_cores
        if use_STAR:
            nc = 1
        parallel.pshell(commands, max_procs=nc, path=out_path)

        # Make sure the analysis completed successfully.
        x1 = [x.log_filename for x in jobs]
        x2 = [x.results_filename for x in jobs]
        x = x1 + x2
        filelib.assert_exists_nz_many(x)
        
        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "rsem"
