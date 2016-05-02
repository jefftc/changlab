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
    # Importing pysam is tricky.  The problem is RSeQC comes with a
    # really old implementation of pysam that does not contain the
    # AlignmentFile object that we need.  Thus, we need to make sure
    # that any other implementation in sys.path is loaded first.
    import sys
    sys_path_old = sys.path[:]
    sys.path = [x for x in sys.path if x.find("RSeQC") < 0]
    import pysam
    sys.path = sys_path_old

    outhandle = open(out_filename, 'w')
    handle = pysam.AlignmentFile(in_filename, "rb")
    
    TAGS = ["MD", "NM", "NH"]
    REQUIRED = ["MD", "NM"]

    header = "query name", "chr", "pos", "CIGAR", "MD", "NM", "NH", \
             "seqlen", "perc_match"
    print >>outhandle, "\t".join(header)
    for i, align in enumerate(handle):
        tag_dict = dict(align.tags)

        # Sometimes reference_id can be -1.  Ignore it.
        if align.reference_id < 0:
            continue

        ref_name = handle.getrname(align.reference_id)
        seqlen = len(align.query_sequence)
        # NH not given by some aligners (e.g. BWA, Bowtie).
        NH = tag_dict.get("NH", "")

        if align.cigarstring is None:
            # No CIGAR string.  Might be unaligned.
            align_pos = ""
            CIGAR = ""
            edit_distance = ""
            perc_match = ""
            MD = ""
            NM = ""
        else:
            missing = [x for x in TAGS if x not in tag_dict]
            missing = [x for x in missing if x in REQUIRED]
            if len(missing) == 1:
                assert not missing, "Missing [%d]: %s tag" % (i, missing[0])
            assert not missing, "Missing tags [%d]: %s" % (
                i, ", ".join(missing))
            #assert "MD" in tag_dict, "Missing: MD tag"
            #assert "NM" in tag_dict, "Missing: NM tag"
            #assert "NH" in tag_dict, "Missing: NH tag"
            align_pos = align.pos
            CIGAR = align.cigarstring
            edit_distance = int(tag_dict["NM"])
            perc_match = 1-float(edit_distance)/seqlen
            MD = tag_dict["MD"]
            NM = tag_dict["NM"]
        x = align.query_name, ref_name, align_pos, CIGAR, MD, NM, NH, \
            seqlen, perc_match
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()
