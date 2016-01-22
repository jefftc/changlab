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

        in_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert in_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in in_filenames:
            p, f = os.path.split(in_filename)
            f, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.intervals" % f)
            x = in_filename, out_filename
            jobs.append(x)
        
        known_sites = []
        x1 = module_utils.get_user_option(
            user_options, "realign_known_sites1", check_file=True)
        x2 = module_utils.get_user_option(
            user_options, "realign_known_sites2", check_file=True)
        x3 = module_utils.get_user_option(
            user_options, "realign_known_sites3", check_file=True)
        x = [x1, x2, x3]
        x = [x for x in x if x]
        known_sites = x
        assert known_sites

        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar -nt 4 \
        #   -T RealignerTargetCreator -R ../genome.idx/erdman.fa -I $i -o $j "

        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, out_filename = x
            x = [("knownSites", x) for x in known_sites]
            x = alignlib.make_GATK_command(
                nt=4, T="RealignerTargetCreator", R=ref.fasta_file_full,
                I=in_filename, o=out_filename, _UNHASHABLE=x)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

    
    def name_outfile(self, antecedents, user_options):
        return "targets"


