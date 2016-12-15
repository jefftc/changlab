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

        bam_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # java -jar picard.jar CollectInsertSizeMetrics \
        #   I=input.bam \
        #   O=insert_size_metrics.txt \
        #   H=insert_size_histogram.pdf \
        #   M=0.5
        opj = os.path.join
        jobs = []   # list of filelib.GenericObject
        for bam_filename in bam_filenames:
            # <in_path>/<sample>.bam
            in_path, sample, ext = mlib.splitpath(bam_filename)
            assert ext == ".bam"

            out_filename = opj(out_path, "%s.insert_size_metrics.txt" % sample)
            hist_filename = opj(
                out_path, "%s.insert_size_histogram.pdf" % sample)
            log_filename = opj(out_path, "%s.log" % sample)
            
            x = filelib.GenericObject(
                sample=sample,
                bam_filename=bam_filename,
                out_filename=out_filename,
                hist_filename=hist_filename,
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
                "-jar", sq(picard_jar), "CollectInsertSizeMetrics",
                "I=%s" % sq(j.bam_filename),
                "O=%s" % sq(j.out_filename),
                "H=%s" % sq(j.hist_filename),
                "M=0.5",
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
        _summarize_insert_size_metrics(jobs, outfile)
        filelib.assert_exists_nz(outfile)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "insert_size"


def _summarize_insert_size_metrics(jobs, outfile):
    # Columns:
    # MEDIAN_INSERT_SIZE
    # MEDIAN_ABSOLUTE_DEVIATION
    # MIN_INSERT_SIZE
    # MAX_INSERT_SIZE
    # MEAN_INSERT_SIZE
    # STANDARD_DEVIATION
    # READ_PAIRS
    # PAIR_ORIENTATION
    # WIDTH_OF_10_PERCENT
    # WIDTH_OF_20_PERCENT
    # WIDTH_OF_30_PERCENT
    # WIDTH_OF_40_PERCENT
    # WIDTH_OF_50_PERCENT
    # WIDTH_OF_60_PERCENT
    # WIDTH_OF_70_PERCENT
    # WIDTH_OF_80_PERCENT
    # WIDTH_OF_90_PERCENT
    # WIDTH_OF_99_PERCENT
    # SAMPLE
    # LIBRARY
    # READ_GROUP

    header = None
    sample2header2value = {}  # sample -> header -> value
    for j in jobs:
        lines = open(j.out_filename).readlines()
        # Find the line that contains this information.
        found = False
        for i, line in enumerate(lines):
            if line.startswith("MEDIAN_INSERT_SIZE"):
                found = True
                break
        assert found, "Missing MEDIAN_INSERT_SIZE: %s" % j.out_filename

        assert i+1 < len(lines)
        x1 = lines[i]
        x2 = lines[i+1]
        x1 = x1.split()
        x2 = x2.split()
        # May be missing values for SAMPLE, LIBRARY, READ_GROUP.
        if len(x2) == len(x1)-3:
            x2 = x2 + ["", "", ""]
        assert len(x1) == len(x2)
        if header is None:
            header = x1
        assert len(x1) == len(header)

        header2value = {}
        for h, v in zip(x1, x2):
            assert h not in header2value
            header2value[h] = v
        sample2header2value[j.sample] = header2value

    handle = open(outfile, 'w')
    all_header = ["Sample"] + header
    print >>handle, "\t".join(all_header)
    for j in jobs:
        header2value = sample2header2value[j.sample]
        values = [header2value[x] for x in header]
        values = [j.sample] + values
        assert len(values) == len(all_header)
        print >>handle, "\t".join(values)
    handle.close()
