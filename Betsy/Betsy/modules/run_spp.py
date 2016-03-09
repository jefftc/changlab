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
            user_options, "control_sample", not_empty=True)

        # Set the experiment name.
        name1 = hashlib.hash_var(treat_sample)
        name2 = hashlib.hash_var(control_sample)
        experiment_name = "%s_vs_%s" % (name1, name2)

        # Make sure the samples exist.
        samples = [x[1] for x in sample_groups]
        assert treat_sample in samples, "Unknown sample: %s" % treat_sample
        assert control_sample in samples, "Unknown sample: %s" % control_sample

        # Find the BAM files.
        treat_filename = run_MACS14.find_bam_file(
            bam_path, treat_sample, sample_groups)
        control_filename = run_MACS14.find_bam_file(
            bam_path, control_sample, sample_groups)
        assert treat_filename, "Missing bam file for %s" % treat_sample
        assert control_filename, "Missing bam file for %s" % control_sample

        cmd = make_pyspp_command(
            treat_filename, control_filename, out_path, num_procs=num_cores)
        log_file = "%s.log" % experiment_name
        cmd = "%s >& %s" % (cmd, log_file)
        parallel.sshell(cmd, path=out_path)

        files = [
            "binding.positions.txt",
            #"broadPeak",
            "crosscorrelation.pdf",
            "density.wig",
            "enrichment.estimates.wig",
            "enrichment.wig",
            #"narrowPeak",   # might be empty if no peaks found
            log_file,
            ]
        filenames = [os.path.join(out_path, x) for x in files]
        filelib.assert_exists_nz_many(filenames)
        
    def name_outfile(self, antecedents, user_options):
        return "spp"


def make_pyspp_command(treat_filename, control_filename, outpath,
                       num_procs=None):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert num_procs is None or (num_procs >= 1 and num_procs < 256)

    # pyspp.py [-j NUM_PROCS] [--fdr_cutoff FDR_CUTOFF]
    #   <treatment_bam> <control_bam> <outpath>
    pyspp = filelib.which_assert(config.pyspp)

    assert os.path.exists(treat_filename)
    assert os.path.exists(control_filename)
    
    sq = parallel.quote
    cmd = [
        sq(pyspp),
        ]
    if num_procs:
        cmd += ["-j", str(num_procs)]
    cmd += [
        sq(treat_filename),
        sq(control_filename),
        sq(outpath),
        ]
    return " ".join(cmd)
