from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils

        bam_filenames = module_utils.find_bam_files(in_data.identifier)
        assert bam_filenames, "No .bam files."
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.matches.txt" % s)
            x = in_filename, out_filename
            jobs.append(x)

        jobs2 = []  # list of (function, args, keywds)
        for x in jobs:
            in_filename, out_filename = x
            x = summarize_bam_file, (in_filename, out_filename), None
            jobs2.append(x)

        parallel.pyfun(jobs2, num_procs=num_cores, DELAY=0.1)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "alignment_matches"


def summarize_bam_file(in_filename, out_filename):
    import pysam

    outhandle = open(out_filename, 'w')

    handle = pysam.AlignmentFile(in_filename, "rb")
    
    header = "query name", "chr", "pos", "CIGAR", "MD", "NM", "NH", \
             "seqlen", "perc_match"
    print >>outhandle, "\t".join(header)
    for i, align in enumerate(handle):
        tag_dict = dict(align.tags)
        assert "MD" in tag_dict, "Missing: MD tag"
        assert "NM" in tag_dict, "Missing: NM tag"
        assert "NH" in tag_dict, "Missing: NH tag"
        ref_name = handle.getrname(align.reference_id)
        seqlen = len(align.query_sequence)
        edit_distance = int(tag_dict["NM"])
        perc_match = 1-float(edit_distance)/seqlen
        x = align.query_name, ref_name, align.pos, align.cigarstring, \
            tag_dict["MD"], tag_dict["NM"], tag_dict["NH"], seqlen, perc_match
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
