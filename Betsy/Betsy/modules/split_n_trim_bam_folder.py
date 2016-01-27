from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

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

        # Make a list of commands.
        sq = shell.quote
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filename = x
            x = alignlib.make_GATK_command(
                T="SplitNCigarReads", R=sq(ref.fasta_file_full),
                I=sq(in_filename), o=sq(out_filename),
                rf="ReassignOneMappingQuality", RMQF=255, RMQT=60,
                U="ALLOW_N_CIGAR_READS")
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        return "split_n_trim.bam"


