from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import parselib
        from genomicode import filelib

        bam_path = in_data.identifier
        assert os.path.exists(bam_path)
        assert os.path.isdir(bam_path)

        # Find all the BAM files.
        bam_filenames = filelib.list_files_in_path(
            bam_path, endswith=".bam", case_insensitive=True)

        results = []  # list of (sample, total, mapped, unmapped, perc mapped)
        for filename in bam_filenames:
            p, f = os.path.split(filename)
            f, e = os.path.splitext(f)
            sample = f

            x = count_mapped_reads(filename)
            mapped, unmapped = x
            total = mapped + unmapped
            perc_mapped = float(mapped)/total*100
            x = sample, total, mapped, unmapped, perc_mapped
            results.append(x)

        # Make table where the rows are the samples and the columns
        # are the statistics.
        table = []
        header = "Sample", "Aligned Reads", "Total Reads", "Perc Aligned"
        table.append(header)
        for x in results:
            sample, total, mapped, unmapped, perc_mapped = x
            
            x1 = parselib.pretty_int(mapped)
            x2 = parselib.pretty_int(total)
            x3 = "%.2f%%" % perc_mapped
            x = sample, x1, x2, x3
            table.append(x)
        
        # Write out the table as text file.
        TXT_FILE = "summary.txt"
        handle = open(TXT_FILE, 'w')
        for x in table:
            print >>handle, "\t".join(x)
        handle.close()

        # XXX txt2xls config
        os.system("txt2xls -b %s > %s" % (TXT_FILE, outfile))
            
        
    def name_outfile(self, antecedents, user_options):
        return "aligned_reads.xls"


def count_mapped_reads(bam_filename):
    # Return a tuple of (mapped, unmapped).
    # bam file must be indexed.
    import os
    import StringIO
    from genomicode import config
    from genomicode import shell
    from genomicode import filelib

    assert os.path.exists(bam_filename)
    
    sq = shell.quote
    samtools = filelib.which_assert(config.samtools)
    cmd = [
        sq(samtools),
        "idxstats",
        sq(bam_filename),
        ]
    # How to check if this is broken?
    x = shell.single(cmd)
    handle = StringIO.StringIO(x)
    total_mapped = total_unmapped = 0
    for x in filelib.read_cols(handle):
        seqname, length, mapped, unmapped = x
        total_mapped += int(mapped)
        total_unmapped += int(unmapped)
    return total_mapped, total_unmapped
    
