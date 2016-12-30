from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        bam_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        metadata = {}
        metadata["tool"] = "samtools %s" % alignlib.get_samtools_version()

        jobs = []
        for bam_filename in bam_filenames:
            x = count_duplicates, (bam_filename,), {}
            jobs.append(x)
        results = parallel.pyfun(jobs, num_procs=num_cores)
        metadata["num_cores"] = num_cores
        assert len(results) == len(bam_filenames)

        handle = open(outfile, 'w')
        header = "Sample", "Duplicated Reads", "Total Reads", "% Duplicated"
        print >>handle, "\t".join(header)
        for i in range(len(bam_filenames)):
            x, sample, x = mlib.splitpath(bam_filenames[i])
            total_reads, dup_reads = results[i]
            perc_dup = float(dup_reads) / total_reads * 100
            perc_dup = "%.2f" % perc_dup
            x = sample, dup_reads, total_reads, perc_dup
            print >>handle, "\t".join(map(str, x))
       
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "duplication.txt"


def count_duplicates(bam_filename):
    # Return a tuple of (total_reads, duplicated_reads).
    import subprocess
    from genomicode import samtools
    from Betsy import module_utils as mlib

    samtools_bin = mlib.get_config("samtools", which_assert_file=True)
    cmd = [
        samtools_bin,
        "view",
        bam_filename,
        ]
    total = num_dup = 0
    p = subprocess.Popen(
        cmd, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    w, r = p.stdin, p.stdout
    w.close()
    for align in samtools.parse_sam(r):
        if align.flag & 0x400:
            num_dup += 1
        total += 1
    return total, num_dup
    
                                
