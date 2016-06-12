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
        genome_size = module_utils.get_user_option(
            user_options, "macs_genome", not_empty=True)
        x = module_utils.get_user_option(
            user_options, "broad_peaks", allowed_values=["no", "yes"])
        broad_peaks = (x == "yes")
        x = module_utils.get_user_option(
            user_options, "macs_paired", allowed_values=["no", "yes"])
        is_paired = (x == "yes")

        # Set the name.
        name = hashlib.hash_var(treat_sample)
        if control_sample:
            x = hashlib.hash_var(control_sample)
            name = "%s_vs_%s" % (treat_sample, x)

        # Make sure the samples exist.
        samples = [x[1] for x in sample_groups]
        assert treat_sample in samples, "Unknown sample: %s" % treat_sample
        if control_sample:
            assert control_sample in samples, \
                   "Unknown sample: %s" % control_sample

        # Find the BAM files.
        treat_filename = run_MACS14.find_bam_file(
            bam_path, treat_sample, sample_groups)
        assert treat_filename, "Missing bam file for %s" % treat_sample
        control_filename = None
        if control_sample:
            control_filename = run_MACS14.find_bam_file(
                bam_path, control_sample, sample_groups)
            assert control_filename, "Missing bam file for %s" % control_sample

        cmd = make_macs2_command(
            treat_filename, control_filename=control_filename, 
            genome_size=genome_size, save_bedgraph_file=True, name=name,
            normalize_read_counts=True, paired=is_paired,
            broad_peak_calling=broad_peaks)
        parallel.sshell(cmd, path=out_path)

        files = [
            "%s_peaks.xls" % name,
            ]
        filenames = [os.path.join(out_path, x) for x in files]
        filelib.assert_exists_nz_many(filenames)
        
    def name_outfile(self, antecedents, user_options):
        return "macs21"


def make_macs2_command(
    treat_filename, control_filename=None, genome_size=None, name=None,
    save_bedgraph_file=False, broad_peak_calling=False, 
    normalize_read_counts=False, paired=False):
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert genome_size in ["hs", "mm", "ce", "dm"]

    # Regular peak calling:
    # macs2 callpeak -t sample.bam -c control.bam \
    #   -f [BAM,BAMPE] -g hs -n T_53BP1 -B -q 0.01
    # 
    # Broad peak calling:
    # macs2 callpeak --broad -t sample.bam -c control.bam \
    #   -f [BAM,BAMPE] -g hs -n T_53BP1 --broad-cutoff 0.1
    # 
    # -n  name.  For saving output files.
    # -w  Save extended fragment pileup at every WIGEXTEND bp in wiggle
    #     file.
    # -B  Save extended fragment pileup at every bp in a bedGraph file.
    #     Much smaller than wiggle file.
    # --broad-cutoff  q-value for merging broad regions.
    # --SPMR          Normalize coverage plot by millions of reads.
    macs2 = filelib.which_assert(config.macs2)
    
    sq = parallel.quote
    cmd = [
        sq(macs2),
        "callpeak",
        ]
    if broad_peak_calling:
        cmd += ["--broad"]
    if normalize_read_counts:
        cmd += ["--SPMR"]
    cmd += ["-t", sq(treat_filename)]
    if control_filename:
        cmd += ["-c", sq(control_filename),]
    format_ = "BAM"
    if paired:
        format_ = "BAMPE"
    cmd += [
        "-f", format_,
        "-g", genome_size,
        ]
    if name:
        cmd.extend(["-n", sq(name)])
    if save_bedgraph_file:
        cmd.append("-B")
    return " ".join(cmd)
