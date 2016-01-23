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
        
        bam_node, ref_node, target_node = antecedents

        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        target_filenames = filelib.list_files_in_path(
            target_node.identifier, endswith=".intervals")
        assert target_filenames, "No .intervals files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        assert len(bam_filenames) == len(target_filenames), \
               "Should have an .intervals file for each bam file."
        sample2bamfilename = {}
        for filename in bam_filenames:
            p, f = os.path.split(filename)
            sample, ext = os.path.splitext(f)
            assert sample not in sample2bamfilename
            sample2bamfilename[sample] = filename
        sample2targetfilename = {}
        for filename in target_filenames:
            p, f = os.path.split(filename)
            sample, ext = os.path.splitext(f)
            assert sample not in sample2targetfilename
            sample2targetfilename[sample] = filename
        assert len(sample2bamfilename) == len(sample2targetfilename)

        missing = [
            x for x in sample2bamfilename if x not in sample2targetfilename]
        assert not missing, "Missing interval files for %d bam files." % \
               len(missing)

        # list of (bam_filename, target_filename, log_filename, out_filename)
        jobs = []
        for sample in sample2bamfilename:
            bam_filename = sample2bamfilename[sample]
            target_filename = sample2targetfilename[sample]
            
            p, f = os.path.split(bam_filename)
            sample, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.bam" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = bam_filename, target_filename, log_filename, out_filename
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

        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar \
        #   -T IndelRealigner -R <ref.fa> \
        #   -I <bam_file> -targetIntervals <target_file> -o <bam_file>
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            bam_filename, target_filename, log_filename, out_filename = x
            x = [("known", x) for x in known_sites]
            x = alignlib.make_GATK_command(
                T="IndelRealigner", R=ref.fasta_file_full,
                I=bam_filename, targetIntervals=target_filename,
                o=out_filename, _UNHASHABLE=x)
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "realigned.bam"


