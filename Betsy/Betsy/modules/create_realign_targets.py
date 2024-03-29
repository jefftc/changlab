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
        
        bam_node, ref_node = antecedents
        in_filenames = mlib.find_bam_files(bam_node.identifier)
        assert in_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        jobs = []  # list of (in_filename, log_filename, out_filename)
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            f, ext = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % f)
            out_filename = os.path.join(out_path, "%s.intervals" % f)
            x = in_filename, log_filename, out_filename
            jobs.append(x)

        filter_reads_with_N_cigar = mlib.get_user_option(
            user_options, "filter_reads_with_N_cigar",
            allowed_values=["no", "yes"])

        known_sites = []
        x1 = mlib.get_user_option(
            user_options, "realign_known_sites1", check_file=True)
        x2 = mlib.get_user_option(
            user_options, "realign_known_sites2", check_file=True)
        x3 = mlib.get_user_option(
            user_options, "realign_known_sites3", check_file=True)
        x = [x1, x2, x3]
        x = [x for x in x if x]
        known_sites = x
        assert known_sites

        # I/O bound, so not likely to get a big speedup with nt.

        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar -nt 4
        #   -T RealignerTargetCreator -R ../genome.idx/erdman.fa -I $i -o $j
        #   --known <known_vcf_file>

        # RealignerTargetCreator takes ~10Gb per process.  Each thread
        # takes the full amount of memory.
        nc = mlib.calc_max_procs_from_ram(12, upper_max=num_cores)

        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filename = x
            
            n = max(1, nc/len(jobs))
            x = [("-known", x) for x in known_sites]
            if filter_reads_with_N_cigar == "yes":
                x.append(("-filter_reads_with_N_cigar", None))
            x = alignlib.make_GATK_command(
                nt=n, T="RealignerTargetCreator", R=ref.fasta_file_full,
                I=in_filename, o=out_filename,
                _UNHASHABLE=x)
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=nc)
        metadata["num_procs"] = nc
        metadata["commands"] = commands

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)
        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "targets"


