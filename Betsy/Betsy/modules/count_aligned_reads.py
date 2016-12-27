from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import parselib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        MAX_CORES = 4  # I/O intensive.

        fastq_node, sample_node, bam_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        sample2fastq = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier, as_dict=True)

        metadata = {}

        jobs = []  # list of (sample, bam_file, fastq_file)
        for filename in bam_filenames:
            path, sample, ext = mlib.splitpath(filename)
            assert sample in sample2fastq, "Missing fastq: %s" % sample
            fastq1, fastq2 = sample2fastq[sample]
            x = sample, filename, fastq1
            jobs.append(x)

        funcalls = []
        for x in jobs:
            sample, bam_filename, fastq_filename = x
            # Count the number of reads.
            x1 = count_reads, (fastq_filename,), {}
            # Count the number of alignments.
            x2 = count_alignments, (bam_filename,), {}
            funcalls.append(x1)
            funcalls.append(x2)
        assert len(funcalls) == len(jobs)*2

        nc = min(num_cores, MAX_CORES)
        results = parallel.pyfun(funcalls, num_procs=nc)
        metadata["num_cores"] = nc

        # list of (sample, aligns, aligned_reads, total_reads, perc_aligned).
        results2 = []
        for i, x in enumerate(jobs):
            sample, bam_filename, fastq_filename = x
            x1 = results[i*2]
            x2 = results[i*2+1]
            total_reads = x1
            aligned_reads, alignments = x2
            perc_aligned = float(aligned_reads)/total_reads
            x = sample, alignments, aligned_reads, total_reads, perc_aligned
            results2.append(x)
        results = results2

        # sort by sample name
        results.sort()

        # Make table where the rows are the samples and the columns
        # are the statistics.
        table = []
        header = ("Sample", "Alignments", "Aligned Reads", "Total Reads",
                  "Perc Aligned")
        table.append(header)
        for x in results:
            sample, alignments, aligned_reads, total_reads, perc_aligned = x

            x1 = parselib.pretty_int(alignments)
            x2 = parselib.pretty_int(aligned_reads)
            x3 = parselib.pretty_int(total_reads)
            x4 = "%.2f%%" % (perc_aligned*100)
            x = sample, x1, x2, x3, x4
            assert len(x) == len(header)
            table.append(x)
        
        # Write out the table as text file.
        TXT_FILE = "summary.txt"
        handle = open(TXT_FILE, 'w')
        for x in table:
            print >>handle, "\t".join(x)
        handle.close()

        txt2xls = mlib.findbin("txt2xls", quote=True)
        parallel.sshell("%s -b %s > %s" % (txt2xls, TXT_FILE, outfile))
        return metadata
            
        
    def name_outfile(self, antecedents, user_options):
        return "aligned_reads.xls"


def count_reads(fastq_filename):
    # Requires an uncompressed fastq file.
    from genomicode import filelib
    from genomicode import parallel

    sq = parallel.quote

    # Make sure it's a fastq file.
    # @M03807:17:000000000-AHGYH:1:1101:20554:1508 1:N:0:16
    # CTTTACACCCAGTGGAGAAGCTCCCAACCAAGCTCTCTTGAGGATCTTGAAGGAAACTGA
    # +
    # <BCC@FAFEC8,C<8968<@EEEFFCCFEC@EDEFGGGGA,@,@EFGGF9,,88,@FFA<
    handle = filelib.openfh(fastq_filename)
    x = [handle.readline() for i in range(4)]
    x = [x.strip() for x in x]
    x = [x for x in x]
    assert len(x) == 4
    assert len(x[1]) == len(x[3])
    assert x[2] == "+"
    
    wc_out = parallel.sshell("wc -l %s" % sq(fastq_filename))
    # velocitron:biocore$ wc -l test01.txt
    # 22278 test01.txt
    # 0 test 1.txt
    x = wc_out.strip().split()
    assert len(x) >= 2, "Unknown format from wc -l\n" % wc_out
    num_lines, filename = x[0], " ".join(x[1:])
    num_lines = int(num_lines)
    num_reads = num_lines / 4
    return num_reads


def count_alignments(bam_filename):
    import subprocess
    from Betsy import module_utils as mlib

    samtools = mlib.findbin("samtools")
    x = [samtools, "view", bam_filename]
    p = subprocess.Popen(
        x, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    r = p.stdout

    alignments = 0
    aligned_reads = {}
    for line in r:
        # M03807:17:000000000-AHGYH:1:2108:11122:14861 99 1 14172 0 12S128...
        # ST-J00106:110:H5NY5BBXX:4:1101:1702:1209 2:N:0:NTCACG   141     *
        x = line.split("\t")
        assert len(x) >= 11, "SAM format"
        qname, flag = x[:2]
        flag = int(flag)
        # 2  mapped in proper pair
        # 4  query is unmapped
        # 8  mate is unmapped
        #is_aligned = flag & 0x02
        is_aligned = not (flag & 0x04)
        if is_aligned:
            alignments += 1
            aligned_reads[qname] = 1
    aligned_reads = len(aligned_reads)
    return aligned_reads, alignments


def count_mapped_reads(bam_filename):
    # Return a tuple of (mapped, unmapped).
    # bam file must be indexed.
    from genomicode import alignlib

    total_mapped = total_unmapped = 0
    for x in alignlib.call_samtools_idxstats(bam_filename):
        seqname, length, mapped, unmapped = x
        total_mapped += int(mapped)
        total_unmapped += int(unmapped)
    return total_mapped, total_unmapped
    
