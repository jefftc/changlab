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
        from Betsy import module_utils

        fastq_node, sample_node, align_node = antecedents
        fastq_data = module_utils.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        assert fastq_data, "I could not find any FASTQ files."
        align_filenames = filelib.list_files_in_path(
            align_node.identifier, endswith=".matches.txt")
        assert align_filenames, "No .matches.txt files."
        align_filenames.sort()

        assert len(fastq_data) == len(align_filenames), \
               "Mismatch: num samples %d %d" % (
            len(fastq_data), len(align_filenames))


        sample2fastqdata = {}
        for x in fastq_data:
            sample, f1, f2 = x
            sample2fastqdata[sample] = x
        
        # list of (sample, align_filename, fastq_filename1, fastq_filename2)
        jobs = []
        for in_filename in align_filenames:
            p, f = os.path.split(in_filename)
            # <sample>.matches.txt
            ext = ".matches.txt"
            assert f.endswith(ext)
            sample = f[:-len(ext)]
            assert sample in sample2fastqdata, "Missing FASTQ: %s" % sample
            x, fastq_filename1, fastq_filename2 = sample2fastqdata[sample]
            x = sample, in_filename, fastq_filename1, fastq_filename2
            jobs.append(x)

        jobs2 = []  # list of (function, args, keywds)
        for x in jobs:
            sample, align_filename, fastq_file1, fastq_file2 = x
            x = align_filename, fastq_file1, fastq_file2
            x = summarize_matches_file, x, None
            jobs2.append(x)

        results = parallel.pyfun(jobs2, num_procs=num_cores, DELAY=0.1)
        assert len(results) == len(jobs2)

        # Put together the results in a table.
        handle = open(out_filename, 'w')
        header = "sample", "perfect", "total", "perc"
        print >>handle, "\t".join(header)
        for x in zip(jobs, results):
            x, d = x
            sample, in_filename, fastq_filename1, fastq_filename2 = x
            x = sample, d["perfect_alignments"], d["total_alignments"], \
                d["perc_perfect"]
            print >>handle, "\t".join(map(str, x))
        handle.close()
    
    def name_outfile(self, antecedents, user_options):
        return "perfect_alignments.txt"


def summarize_matches_file(filename, fastq_file1, fastq_file2):
    # Return dictionary with keys:
    # total_alignments       int
    # perfect_alignments     int
    # perc_perfect           float (0.0-1.0)
    from genomicode import filelib
    from genomicode import genomelib

    # Get the list of all possible read names from fastq_file1.
    # Title from fastq file:
    #   @ST-J00106:107:H5NK2BBXX:1:1101:1438:1173 1:N:0:NAGATC
    # From alignment file:
    #   ST-J00106:107:H5NK2BBXX:1:2218:22079:11653
    all_aligns = {}
    for x in genomelib.read_fastq(fastq_file1):
        title, sequence, quality = x
        x = title
        if x.startswith("@"):
            x = x[1:]
        x = x.split()[0]  # alignment file only contains the first part.
        all_aligns[x] = 1
    
    perfect_aligns = {}
    for d in filelib.read_row(filename, header=1):
        assert d.query_name in all_aligns
        if int(d.NM) == 0:
            perfect_aligns[d.query_name] = 1

    perfect = len(perfect_aligns)
    total = len(all_aligns)
    results = {
        "perfect_alignments" : perfect,
        "total_alignments" : total,
        "perc_perfect" : float(perfect) / total,
        }
    return results
