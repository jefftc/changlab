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

        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # java -jar picard.jar CollectAlignmentSummaryMetrics \
        #   R=reference_sequence.fasta \
        #   I=input.bam \
        #   O=output.txt
        opj = os.path.join
        jobs = []   # list of filelib.GenericObject
        for bam_filename in bam_filenames:
            # <in_path>/<sample>.bam
            in_path, sample, ext = mlib.splitpath(bam_filename)
            assert ext == ".bam"
            out_filename = opj(out_path, "%s.alignment_metrics.txt" % sample)
            log_filename = opj(out_path, "%s.log" % sample)
            x = filelib.GenericObject(
                sample=sample,
                bam_filename=bam_filename,
                out_filename=out_filename,
                log_filename=log_filename)
            jobs.append(x)

        # Make the commands to run picard.
        picard_jar = alignlib.find_picard_jar("picard")
        sq = parallel.quote
        commands = []
        for j in jobs:
            # Should have better way of getting java path.
            cmd = [
                "java",
                "-Xmx10g",
                "-jar", sq(picard_jar), "CollectAlignmentSummaryMetrics",
                "I=%s" % sq(j.bam_filename),
                "R=%s" % sq(ref.fasta_file_full),
                "O=%s" % sq(j.out_filename),
                ]
            cmd = " ".join(cmd)
            cmd = "%s >& %s" % (cmd, sq(j.log_filename))
            commands.append(cmd)

        metadata["commands"] = commands
        parallel.pshell(commands, max_procs=num_cores)
        x = [x.out_filename for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Summarize the insert size files.
        outfile = opj(out_path, "summary.txt")
        _summarize_alignment_summary_metrics(jobs, outfile)
        filelib.assert_exists_nz(outfile)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "alignment"


def _summarize_alignment_summary_metrics(jobs, outfile):
    # Columns:
    # CATEGORY       FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR  (what about SE?)
    # TOTAL_READS
    # PF_READS
    # PCT_PF_READS
    # PF_NOISE_READS
    # PF_READS_ALIGNED
    # PCT_PF_READS_ALIGNED
    # PF_ALIGNED_BASES
    # PF_HQ_ALIGNED_READS
    # PF_HQ_ALIGNED_BASES
    # PF_HQ_ALIGNED_Q20_BASES
    # PF_HQ_MEDIAN_MISMATCHES
    # PF_MISMATCH_RATE
    # PF_HQ_ERROR_RATE
    # PF_INDEL_RATE
    # MEAN_READ_LENGTH
    # READS_ALIGNED_IN_PAIRS
    # PCT_READS_ALIGNED_IN_PAIRS
    # BAD_CYCLES
    # STRAND_BALANCE
    # PCT_CHIMERAS
    # PCT_ADAPTER
    # SAMPLE
    # LIBRARY
    # READ_GROUP
    handle = open(outfile, 'w')
    header = None
    header_printed = False
    for j in jobs:
        for line in open(j.out_filename):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            cols = line.split()
            if not header:
                 header = cols
            if line.startswith("CATEGORY"):
                assert header == cols
                if not header_printed:
                    header_printed = True
                    cols = ["Sample"] + cols
                    print >>handle, "\t".join(cols)
            else:
                # May be missing values for SAMPLE, LIBRARY, READ_GROUP.
                if len(cols) == len(header)-3:
                    cols = cols + ["", "", ""]
                assert len(cols) == len(header)
                cols = [j.sample] + cols
                print >>handle, "\t".join(cols)
