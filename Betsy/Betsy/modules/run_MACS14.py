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
        from genomicode import config
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
        shiftsize = module_utils.get_user_option(
            user_options, "macs_shiftsize")
        if shiftsize:
            shiftsize = int(shiftsize)

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
            treat_filename, control_filename, name=name,
            genome_size=genome_size, shiftsize=shiftsize,
            save_bedgraph_file=True)
        shell.single(cmd, path=out_path)

        # Run Rscript on the model, if one was generated.
        model_file = os.path.join(out_path, "%s_model.r" % name)
        if os.path.exists(model_file):
            Rscript = filelib.which_assert(config.Rscript)
            cmd = [shell.quote(Rscript), model_file]
            shell.single(cmd, path=out_path)

        files = [
            "%s_peaks.xls" % name,
            "%s_summits.bed" % name,
            ]
        filelib.assert_exists_nz_many(files)

        
    def name_outfile(self, antecedents, user_options):
        return "macs14"


def make_macs14_command(
    treat_filename, control_filename=None, genome_size=None, name=None,
    shiftsize=None, save_wiggle_file=False, save_single_wiggle_file=False,
    save_bedgraph_file=False, call_subpeaks=False):
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert genome_size in ["hs", "mm", "ce", "dm"]
    if call_subpeaks:
        save_wiggle_file = True
        save_bedgraph_file = False
    if shiftsize:
        assert shiftsize > 0 and shiftsize < 10000

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
    #
    # If estimated fragment size is too short (e.g. 53), then specify
    # your own fragment size.  shiftsize is 1/2 of fragment size.
    # --nomodel --shiftsize 73 (for fragment size of 146)
    # Often fragment size is 150-200 for ChIP-Seq.
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
    if shiftsize:
        cmd.extend([
            "--nomodel",
            "--shiftsize", str(shiftsize),
            ])
    if save_wiggle_file:
        cmd.append("-w")
    if save_single_wiggle_file:
        cmd.append("-S")
    if save_bedgraph_file:
        cmd.append("-B")
    if call_subpeaks:
        cmd.append("--call_subpeaks")
    return " ".join(cmd)
