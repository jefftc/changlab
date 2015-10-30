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
            user_options, "treatment_sample", not_empty=True)
        control_sample = module_utils.get_user_option(
            user_options, "control_sample")
        genome_size = module_utils.get_user_option(
            user_options, "macs_genome", not_empty=True)
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
        treat_filename = find_bam_file(bam_path, treat_sample, sample_groups)
        assert treat_filename, "Missing bam file for %s" % treat_sample
        control_filename = None
        if control_sample:
            control_filename = find_bam_file(
                bam_path, control_sample, sample_groups)
            assert control_filename, "Missing bam file for %s" % control_sample

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
        filenames = [os.path.join(out_path, x) for x in files]
        filelib.assert_exists_nz_many(filenames)

        
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


def find_bam_file(path, sample, sample_groups):
    # Return the full path of a bam file or None, if not found.
    
    # Try looking for:
    # 1.  <path>/<sample>.bam
    #     This is generated by the alignment pipeline.
    # 2.  <path>/<filename>.bam in sample_groups
    #     This might be given by the user.
    #     The sample group file may contain either FASTA or BAM files.
    import os

    filename = os.path.join(path, "%s.bam" % sample)
    if os.path.exists(filename):
        return filename
    
    # Look in sample_groups.
    # <path>/<filename>
    # <path>/<filestem>.bam

    # <filestem>
    # Have to remove:
    # 1.  COMPRESS_EXTS
    # 2.  Other extensions (.fasta, .fa)
    # 3.  _1, _2
    COMPRESS_EXTS = [".gz", ".bz2", ".xz"]
    
    for x in sample_groups:
        filename, s, pair = x
        if s != sample:
            continue
        if filename.lower().endswith(".bam"):
            x = os.path.join(path, filename)
            if os.path.exists(x):
                return x
        
        # Make filestems.
        filestems = []
        fs = filename
        # No compression extensions.
        for x in COMPRESS_EXTS:
            if fs.lower().endswith(x):
                fs = fs[:len(x)]
        filestems.append(fs)
        # No extensions at all.
        fs, e = os.path.splitext(fs)
        filestems.append(fs)
        # No _1 or _2
        if fs.endswith("_1") or fs.endswith("_2"):
            fs = fs[:-2]
        filestems.append(fs)

        # See if any of the filestems can be made into bam files.
        for fs in filestems:
            x = "%s.bam" % fs
            x = os.path.join(path, x)
            if os.path.exists(x):
                return x
        
    return None
    
    
