from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import shell
        from genomicode import hashlib
        from genomicode import filelib
        from Betsy import module_utils

        bam_node, group_node = antecedents
        bam_path = module_utils.check_inpath(bam_node.identifier)
        sample_groups = module_utils.read_sample_group_file(
            group_node.identifier)

        # Get options.
        treat_sample = module_utils.get_user_option(
            user_options, "treatment_sample", required=True)
        control_sample = module_utils.get_user_option(
            user_options, "control_sample")
        genome_size = module_utils.get_user_option(
            user_options, "macs_genome", required=True)
        x = module_utils.get_user_option(
            user_options, "broad_peaks", default="no",
            allowed_values=["no", "yes"])
        broad_peaks = (x == "yes")

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

        # Assume if one sample is paired, all samples are paired.
        pairs = [x[2] for x in sample_groups]
        pairs = [x for x in x if x != None]
        is_paired = False
        if pairs:
            is_paired = True

        # Find the BAM files.
        # Should be named:
        # <sample>.bam
        treat_file = "%s.bam" % treat_sample
        treat_filename = os.path.join(bam_path, treat_file)
        assert os.path.exists(treat_filename), \
               "File not found: %s" % treat_file

        control_filename = None
        if control_sample:
            control_file = "%s.bam" % control_sample
            control_filename = os.path.join(bam_path, control_file)
            assert os.path.exists(control_filename), \
                   "File not found: %s" % control_file

        cmd = make_macs2_command(
            treat_filename, control_filename, 
            genome_size=genome_size, save_bedgraph_file=True, name=name,
            normalize_read_counts=True, paired=is_paired,
            broad_peak_calling=broad_peaks)
        shell.single(cmd, path=out_path)

        files = [
            "%s_peaks.xls" % name,
            ]
        filelib.assert_exists_nz_many(files)
        
    def name_outfile(self, antecedents, user_options):
        return "macs21"


def make_macs2_command(
    treat_filename, control_filename=None, genome_size=None, name=None,
    save_bedgraph_file=False, broad_peak_calling=False, 
    normalize_read_counts=False, paired=False):
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

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
    
    sq = shell.quote
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
