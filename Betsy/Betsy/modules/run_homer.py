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
        from genomicode import config
        from Betsy import module_utils

        tag_node, group_node = antecedents
        tag_path = module_utils.check_inpath(tag_node.identifier)
        sample_groups = module_utils.read_sample_group_file(
            group_node.identifier)

        # Get options.
        treat_sample = module_utils.get_user_option(
            user_options, "treatment_sample", not_empty=True)
        control_sample = module_utils.get_user_option(
            user_options, "control_sample")

        # Set the experiment name.
        experiment_name = treat_sample
        if control_sample:
            name1 = hashlib.hash_var(treat_sample)
            name2 = hashlib.hash_var(control_sample)
            experiment_name = "%s_vs_%s" % (name1, name2)

        # Make sure the samples exist.
        samples = [x[1] for x in sample_groups]
        assert treat_sample in samples, "Unknown sample: %s" % treat_sample
        assert control_sample in samples, "Unknown sample: %s" % control_sample

        # Find the tag directories.
        treat_path = os.path.join(tag_path, treat_sample)
        assert os.path.exists(treat_path)
        if control_sample:
            control_path = os.path.join(tag_path, control_sample)
            assert os.path.exists(control_path)

        # Get the command.
        homer_path = filelib.which_assert(config.homer_path)
        x = os.path.join(homer_path, "bin", "findPeaks")
        assert filelib.exists_nz(x)
        find_peaks = x

        log_file = "%s.log" % experiment_name
        peak_file = "%s.peaks.txt" % experiment_name

        sq = parallel.quote
        cmd = [
            sq(find_peaks),
            sq(treat_path),
            "-style", "factor",
            ]
        if control_sample:
            cmd += ["-i", control_path]
        cmd = " ".join(cmd)
        cmd = "%s 2> %s 1> %s" % (cmd, log_file, peak_file)
        parallel.sshell(cmd, path=out_path)

        x = os.path.join(out_path, peak_file)
        filelib.assert_exists_nz(x)

        
    def name_outfile(self, antecedents, user_options):
        return "homer"
    
