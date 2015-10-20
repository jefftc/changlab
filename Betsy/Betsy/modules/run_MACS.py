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
        from Betsy import module_utils

        bam_node, group_node = antecedents
        bam_path = bam_node.identifier
        assert os.path.exists(bam_path)
        assert os.path.isdir(bam_path)
        filelib.safe_mkdir(out_path)

        # Figure out which samples to process.
        
        TREATMENT = "treatment_sample"
        assert TREATMENT in user_options, "Missing option: %s" % TREATMENT
        treat_sample = user_options[TREATMENT]
        assert treat_sample, "Empty treatment sample."

        # CONTROL is optional.
        CONTROL = "control_sample"
        control_sample = None
        if CONTROL in user_options:
            control_sample = user_options[CONTROL]
            assert control_sample, "Empty control sample."

        GENOME = "macs_genome"
        assert GENOME in user_options, "Missing option: %s" % GENOME
        genome_size = user_options[GENOME]

        # Make sure the samples exist.
        x = module_utils.read_sample_group_file(group_node.identifier)
        sample_groups = x
        samples = [x[1] for x in sample_groups]
        assert treat_sample in samples, "Unknown sample: %s" % treat_sample
        if control_sample:
            assert control_sample in samples, \
                   "Unknown sample: %s" % control_sample

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

        cmd = make_macs14_command(
            treat_filename, control_filename, 
            genome_size=genome_size, save_bedgraph_file=True)
        shell.single(cmd)

        
    def name_outfile(self, antecedents, user_options):
        return "MACS"


def make_macs14_command(
    treat_filename, control_filename=None, genome_size=None, name=None,
    save_wiggle_file=False, save_single_wiggle_file=False,
    save_bedgraph_file=False, call_subpeaks=False):
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert genome_size in ["hs", "mm", "ce", "dm"]
    if call_subpeaks:
        save_wiggle_file = True
        save_bedgraph_file = False

    #macs14 -t Sample_4_T_53BP1.sorted.bam -c Sample_8_T_input.sorted.bam \
    #   -g hs -n T_53BP1 -B -S --call-subpeaks >& T_53BP1.log
    # -n  name.  For saving output files.
    # -w  Save extended fragment pileup at every WIGEXTEND bp in wiggle
    #     file.
    # -B  Save extended fragment pileup at every bp in a bedGraph file.
    #     Much smaller than wiggle file.
    # -S  A single wiggle file will be saved for treatment and input.
    #     i.e. for whole genome, rather than for each chromosome.
    # --call_subpeaks  Use PeakSplitter algorithm to find subpeaks.
    #                  -w needs to be on, and -B should be off.
    macs14 = filelib.which_assert(config.macs14)
    
    sq = shell.quote
    cmd = [
        sq(macs14),
        "-t", sq(treat_filename),
        ]
    if control_filename:
        cmd += ["-c", sq(control_filename),]
    cmd += [
        "-f", "BAM",
        "-g", genome_size,
        ]
    if name:
        cmd.extend(["-n", sq(name)])
    if save_wiggle_file:
        cmd.append("-w")
    if save_single_wiggle_file:
        cmd.append("-S")
    if save_bedgraph_file:
        cmd.append("-B")
    if call_subpeaks:
        cmd.append("--call_subpeaks")
    return " ".join(cmd)
