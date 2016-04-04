from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        fastq_node, sample_node, align_node = antecedents
        fastq_data = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        assert fastq_data, "I could not find any FASTQ files."
        align_filenames = filelib.list_files_in_path(
            align_node.identifier, endswith=".matches.txt")
        assert align_filenames, "No .matches.txt files."
        align_filenames.sort()
        metadata = {}

        assert len(fastq_data) == len(align_filenames), \
               "Mismatch: num samples %d %d" % (
            len(fastq_data), len(align_filenames))

        num_mismatches = mlib.get_user_option(
            user_options, "num_mismatches", type=int)
        assert num_mismatches >= 0 and num_mismatches < 25
        metadata["num_mismatches"] = num_mismatches

        sample2fastqdata = {}
        for x in fastq_data:
            sample, f1, f2 = x
            sample2fastqdata[sample] = x
        
        # list of (sample, align_filename, summary_filename,
        #   fastq_filename1, fastq_filename2)
        jobs = []
        for in_filename in align_filenames:
            p, f = os.path.split(in_filename)
            # <sample>.matches.txt
            ext = ".matches.txt"
            assert f.endswith(ext)
            sample = f[:-len(ext)]
            assert sample in sample2fastqdata, "Missing FASTQ: %s" % sample
            summary_filename = "%s.summary.txt" % sample
            x, fastq_filename1, fastq_filename2 = sample2fastqdata[sample]
            x = sample, in_filename, summary_filename, \
                fastq_filename1, fastq_filename2
            jobs.append(x)

        jobs2 = []  # list of (function, args, keywds)
        for x in jobs:
            sample, align_filename, summary_filename, \
                    fastq_file1, fastq_file2 = x
            args = align_filename, fastq_file1, fastq_file2, num_mismatches
            keywds = {
                "temp_path" : ".",
                "outfile" : summary_filename,
                }
            x = summarize_matches_file, args, keywds
            jobs2.append(x)

        # Since this can take a lot of memory (depending on the number
        # of reads, can easily take 8 Gb), do just 1 process at a
        # time.  Also, I/O intensive.  Don't do too many at a time.
        MAX_PROCS = 1
        #MAX_PROCS = 4
        nc = min(MAX_PROCS, num_cores)
        results = parallel.pyfun(jobs2, num_procs=nc, DELAY=0.1)
        metadata["num_cores"] = nc
        assert len(results) == len(jobs2)

        # Put together the results in a table.
        handle = open(out_filename, 'w')
        header = "sample", "match", "total", "perc", "perc mismatch"
        print >>handle, "\t".join(header)
        for x in zip(jobs, results):
            x, d = x
            sample, in_filename, summary_filename, \
                    fastq_filename1, fastq_filename2 = x
            perc_mismatch = 1 - d["perc_perfect"]
            x = sample, d["perfect_alignments"], d["total_alignments"], \
                d["perc_perfect"], perc_mismatch
            assert len(x) == len(header)
            print >>handle, "\t".join(map(str, x))
        handle.close()
        return metadata
    
    def name_outfile(self, antecedents, user_options):
        return "matched_alignments.txt"


def summarize_matches_file(filename, fastq_file1, fastq_file2, num_mismatches,
                           temp_path=None, outfile=None):
    # Return dictionary with keys:
    # total_alignments       int
    # perfect_alignments     int
    # perc_perfect           float (0.0-1.0)
    # Will create temporary files in temp_path.  These files could be
    # big (similar in size to the fastq files), so this should
    # ideally be a path with a lot of free space.
    import os
    import tempfile
    import gdbm
    import json
    from genomicode import filelib
    from genomicode import genomelib

    if temp_path is None:
        temp_path = "."
        
    # Get the list of all possible read names from fastq_file1.
    # Title from fastq file:
    #   @ST-J00106:107:H5NK2BBXX:1:1101:1438:1173 1:N:0:NAGATC
    # From alignment file:
    #   ST-J00106:107:H5NK2BBXX:1:2218:22079:11653

    IN_MEMORY = True

    temp_filename = None
    all_aligns = None
    try:
        if IN_MEMORY:
            # This can take a lot of memory.
            all_aligns = {}
            value = 0
        else:
            # To minimize the use of memory, store in a database.
            x, temp_filename = tempfile.mkstemp(dir=temp_path)
            os.close(x)
            all_aligns = gdbm.open(temp_filename, "nf")
            value = "0"

        for x in genomelib.read_fastq(fastq_file1):
            title, sequence, quality = x
            x = title
            if x.startswith("@"):
                x = x[1:]
            x = x.split()[0]  # alignment file only contains the first part.
            all_aligns[x] = value

        # Keep track of the ones that I've seen before to make sure
        # we don't double count.
        #perfect_aligns = {}
        perfect = 0
        for d in filelib.read_row(filename, header=1):
            # This check makes the function very slow.
            assert d.query_name in all_aligns
            if int(d.NM) <= num_mismatches:
                num_seen = int(all_aligns[d.query_name])
                if num_seen == 0:
                    perfect += 1
                num_seen += 1
                if not IN_MEMORY:
                    num_seen = str(num_seen)
                all_aligns[d.query_name] = num_seen
                #perfect += 1
                #perfect_aligns[d.query_name] = 1
        #perfect = len(perfect_aligns)
        total = len(all_aligns)
    finally:
        if all_aligns and not IN_MEMORY:
            all_aligns.close()
        # dbm will create files:
        # <temp_file>
        # <temp_file>.dir
        # <temp_file>.pag
        # Delete them all.
        if temp_filename is not None:
            temp_file = os.path.split(temp_filename)[1]
            x = os.listdir(temp_path)
            x = [x for x in x if x.startswith(temp_file)]
            x = [os.path.join(temp_path, x) for x in x]
            for x in x:
                os.unlink(x)
        #if os.path.exists(temp_filename):
        #    os.unlink(temp_filename)
    
    results = {
        "perfect_alignments" : perfect,
        "total_alignments" : total,
        "perc_perfect" : float(perfect) / total,
        }
    if outfile is not None:
        x = json.dumps(results)
        open(outfile, 'w').write(x)
    return results
