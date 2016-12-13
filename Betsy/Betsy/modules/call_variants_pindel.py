from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        bam_node, ref_node, insert_size_node, alignment_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # ./pindel -f <reference.fa> -i <bam_configuration_file> 
        #   -c <chromosome_name> -o <out_prefix>
        #   -T <num threads>
        # 
        # Creates files:
        # <out_prefix>_D     Deletion
        # <out_prefix>_SI    Short insertion
        # <out_prefix>_LI    Long insertion
        # <out_prefix>_INV   Inversion
        # <out_prefix>_TD    Tandem deletion
        # <out_prefix>_BP    Breakpoint
        # <out_prefix>_RP    ??? read pair???
        # <out_prefix>_CloseEndMapped   Only on end could be mapped.

        # Pindel cannot handle spaces in the BAM filenames (because of
        # the config file).  Symlink the file to a local directory to make
        # sure there are no spaces.
        bam_path = "bam"

        opj = os.path.join
        jobs = []  # list of filelib.GenericObject
        for bam_filename in bam_filenames:
            p, f = os.path.split(bam_filename)
            sample, ext = os.path.splitext(f)
            bai_filename = "%s.bai" % bam_filename
            filelib.assert_exists_nz(bai_filename)
            x = sample.replace(" ", "_")
            local_bam = opj(bam_path, "%s.bam" % x)
            local_bai = opj(bam_path, "%s.bam.bai" % x)
            config_filename = opj(out_path, "%s.config.txt" % sample)
            out_prefix = opj(out_path, sample)
            log_filename = opj(out_path, "%s.log" % sample)
            x = filelib.GenericObject(
                sample=sample,
                bam_filename=bam_filename, bai_filename=bai_filename,
                local_bam=local_bam, local_bai=local_bai,
                config_filename=config_filename,
                out_prefix=out_prefix,
                log_filename=log_filename)
            jobs.append(x)

        filelib.safe_mkdir(bam_path)
        for j in jobs:
            assert " " not in j.local_bam
            filelib.assert_exists_nz(j.bam_filename)
            filelib.assert_exists_nz(j.bai_filename)
            if not os.path.exists(j.local_bam):
                os.symlink(j.bam_filename, j.local_bam)
            if not os.path.exists(j.local_bai):
                os.symlink(j.bai_filename, j.local_bai)
            
        # Read the insert sizes.
        summary_file = opj(insert_size_node.identifier, "summary.txt")
        filelib.assert_exists_nz(summary_file)
        sample2size = _read_insert_sizes(summary_file)
        # Make sure all the samples have inserts.
        for j in jobs:
            assert j.sample in sample2size, \
                   "Missing in insert size file: %s" % j.sample

        # Read the fragment sizes.
        summary_file = opj(alignment_node.identifier, "summary.txt")
        filelib.assert_exists_nz(summary_file)
        sample2readlen = _read_fragment_sizes(summary_file)
        # Make sure all the samples have read lengths.
        for j in jobs:
            assert j.sample in sample2readlen, \
                   "Missing in alignment summary file: %s" % j.sample

        # Make the config file.
        for j in jobs:
            # <insert size> is the whole length to be sequenced, including
            # the length of the pair of reads.  Picard only counts the
            # sequence between the reads.
            size = sample2size[j.sample]
            read_length = sample2readlen[j.sample]
            insert_size = size + read_length*2
            handle = open(j.config_filename, 'w')
            print >>handle, "%s %s %s" % (j.local_bam, insert_size, j.sample)
            handle.close()

        # Make a list of commands.
        pindel = mlib.get_config("pindel", which_assert_file=True)
        sq = parallel.quote
        commands = []
        for j in jobs:
            cmd = [
                sq(pindel),
                "-f", sq(ref.fasta_file_full),
                "-i", sq(j.config_filename),
                "-c", "ALL",
                "-T", 1,
                "-o", sq(j.out_prefix),
                ]
            cmd = " ".join(map(str, cmd))
            cmd = "%s >& %s" % (cmd, j.log_filename)
            commands.append(cmd)
        parallel.pshell(commands, max_procs=num_cores)
        metadata["num_cores"] = num_cores
        metadata["commands"] = commands

        # Make sure the analysis completed successfully.  If not, try
        # to diagnose.
        x = [x.log_filename for x in jobs]
        filelib.assert_exists_nz_many(x)
        x1 = ["%s_D" % x.out_prefix for x in jobs]
        x2 = ["%s_SI" % x.out_prefix for x in jobs]
        x3 = ["%s_LI" % x.out_prefix for x in jobs]
        x4 = ["%s_INV" % x.out_prefix for x in jobs]
        x5 = ["%s_TD" % x.out_prefix for x in jobs]
        x6 = ["%s_BP" % x.out_prefix for x in jobs]
        x = x1 + x2 + x3 + x4 + x5 + x6
        filelib.assert_exists_many(x)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "pindel"


def _read_insert_sizes(filename):
    from genomicode import filelib

    # Should be summary from Picard's insert size summary.
    sample2size = {}
    for d in filelib.read_row(filename, header=1):
        assert hasattr(d, "Sample"), \
               "Missing in summary: Sample"
        assert hasattr(d, "MEAN_INSERT_SIZE"), \
               "Missing in summary: MEAN_INSERT_SIZE"
        x = float(d.MEAN_INSERT_SIZE)
        size = int(round(x))
        assert size > 0 and size < 10000  # checking
        sample2size[d.Sample] = size
    return sample2size


def _read_fragment_sizes(filename):
    from genomicode import filelib

    # Should be summary from Picard's alignment summary.
    sample2readlen = {}
    for d in filelib.read_row(filename, header=1):
        assert hasattr(d, "Sample"), \
               "Missing in summary: Sample"
        assert hasattr(d, "MEAN_READ_LENGTH"), \
               "Missing in summary: MEAN_READ_LENGTH"
        x = float(d.MEAN_READ_LENGTH)
        readlen = int(round(x))
        assert readlen > 0 and readlen < 10000  # checking
        sample2readlen[d.Sample] = readlen
    return sample2readlen
