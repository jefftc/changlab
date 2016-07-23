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
        from Betsy import module_utils as mlib

        MAX_RAM = 64   # maximum amount of ram to use in Gb.

        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        jobs = []  # list of (in_filename, log_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % s)
            out_filename = os.path.join(out_path, f)
            x = in_filename, log_filename, out_filename
            jobs.append(x)
        
        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar
        #   -T SplitNCigarReads -R ../hg19.fa -I $i -o $j
        #   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
        #   -U ALLOW_N_CIGAR_READS

        # Start with 5 Gb RAM.
        commands = make_commands(jobs, ref.fasta_file_full, 5)
        nc = mlib.calc_max_procs_from_ram(5, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["commands"] = commands
        metadata["num_procs"] = nc

        # If any of the analyses didn't finish, try again with more
        # RAM.
        jobs2 = []
        for x in jobs:
            in_filename, log_filename, out_filename = x
            if filelib.exists_nz(out_filename):
                continue
            jobs2.append(x)
        if jobs2:
            commands = make_commands(jobs2, ref.fasta_file_full, MAX_RAM)
            nc = mlib.calc_max_procs_from_ram(MAX_RAM, upper_max=num_cores)
            parallel.pshell(commands, max_procs=nc)
            metadata["commands"] += commands
            
        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "split_n_trim.bam"


def make_commands(jobs, fasta_filename, memory):
    from genomicode import parallel
    from genomicode import alignlib
    
    sq = parallel.quote
    commands = []
    for x in jobs:
        in_filename, log_filename, out_filename = x
        params = {
            "T" : "SplitNCigarReads",
            "R" : sq(fasta_filename),
            "I" : sq(in_filename),
            "o" : sq(out_filename),
            "rf" : "ReassignOneMappingQuality",
            "RMQF" : 255,
            "RMQT" : 60,
            "U" : "ALLOW_N_CIGAR_READS",
            }
        x = alignlib._make_java_command("gatk_jar", params, 1, memory=memory)
        x = "%s >& %s" % (x, sq(log_filename))
        commands.append(x)
    return commands
    
