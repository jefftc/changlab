from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import hashlib
        from genomicode import filelib
        from Betsy import module_utils
        import run_MACS14

        bam_node, group_node = antecedents
        bam_path = module_utils.check_inpath(bam_node.identifier)
        sample_groups = module_utils.read_sample_group_file(
            group_node.identifier)

        # Get options.
        treat_sample = module_utils.get_user_option(
            user_options, "treatment_sample", not_empty=True)
        control_sample = module_utils.get_user_option(
            user_options, "control_sample")
        fragment_length = module_utils.get_user_option(
            user_options, "peakseq_fragment_length", not_empty=True, type=int)
        mappability_file = module_utils.get_user_option(
            user_options, "mappability_file", not_empty=True, check_file=True)
        assert fragment_length > 0 and fragment_length < 1000


        # Set the experiment name.
        name1 = hashlib.hash_var(treat_sample)
        name2 = hashlib.hash_var(control_sample)
        experiment_name = "%s_vs_%s" % (name1, name2)

        # Make sure the samples exist.
        samples = [x[1] for x in sample_groups]
        assert treat_sample in samples, "Unknown sample: %s" % treat_sample
        if control_sample:
            assert control_sample in samples, \
                   "Unknown sample: %s" % control_sample

        # Find the BAM files.
        treat_filename = run_MACS14.find_bam_file(
            bam_path, treat_sample, sample_groups)
        control_filename = run_MACS14.find_bam_file(
            bam_path, control_sample, sample_groups)
        assert treat_filename, "Missing bam file for %s" % treat_sample
        assert control_filename, "Missing bam file for %s" % control_sample

        cmd = make_peakseq_command(
            treat_filename, control_filename, 
            out_path, experiment_name, fragment_length, mappability_file)
        log_file = "%s.log" % experiment_name
        cmd = "%s >& %s" % (cmd, log_file)
        parallel.sshell(cmd, path=out_path)

        files = [
            "config.dat",
            log_file,
            "%s.txt" % experiment_name,
            # Can be length 0, if no peaks found.
            #"%s_narrowPeak.txt" % experiment_name,
            ]
        filenames = [os.path.join(out_path, x) for x in files]
        filelib.assert_exists_nz_many(filenames)
        
    def name_outfile(self, antecedents, user_options):
        return "peakseq"


def make_peakseq_command(
    treat_filename, control_filename, outpath,
    experiment_name, fragment_length, mappability_file):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    # pypeakseq.py --experiment_name EXPERIMENT_NAME
    #   --fragment_length FRAGMENT_LENGTH
    #   <mappability_file> <treatment_bam> <control_bam> <outpath>
    pypeakseq = filelib.which_assert(config.pypeakseq)

    assert os.path.exists(treat_filename)
    assert os.path.exists(control_filename)
    assert os.path.exists(mappability_file)
    assert fragment_length > 0 and fragment_length < 100000
    
    sq = parallel.quote
    cmd = [
        sq(pypeakseq),
        "--experiment_name", experiment_name,
        "--fragment_length", str(fragment_length),
        sq(mappability_file),
        sq(treat_filename),
        sq(control_filename),
        sq(outpath),
        ]
    return " ".join(cmd)
