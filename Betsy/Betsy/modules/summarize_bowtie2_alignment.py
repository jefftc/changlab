from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib
        from genomicode import parselib

        import summarize_bowtie1_alignment

        log_filenames = summarize_bowtie1_alignment._find_output_logs(
            in_data.identifier)
        assert log_filenames

        results = {}  # dict of sample -> dictionary of output
        for filename in log_filenames:
            # <path>/<sample>.log
            path, file_ = os.path.split(filename)
            f, e = os.path.splitext(file_)
            assert e == ".log"
            sample = f
            results[sample] = module_utils.parse_bowtie2_output(filename)
        all_samples = sorted(results)

        is_paired = False
        if "concordant_reads" in results[all_samples[0]]:
            is_paired = True

        # Make table where the rows are the samples and the columns
        # are the statistics.
        table = []
        header = "Sample", "Aligned Reads", "Total Reads", "Perc Aligned"
        if is_paired:
            header = header + ("Concordant Reads",)
        table.append(header)
        for sample in all_samples:
            stats = results[sample]
            total_reads = stats["reads_processed"]
            aligned_reads = stats["aligned_reads"]
            concordant_reads = stats.get("concordant_reads")
            perc_aligned = float(aligned_reads) / total_reads * 100

            x1 = parselib.pretty_int(aligned_reads)
            x2 = parselib.pretty_int(total_reads)
            x3 = "%.2f%%" % perc_aligned
            x = sample, x1, x2, x3
            if is_paired:
                x4 = parselib.pretty_int(concordant_reads)
                x = x + (x4,)
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
        return "bowtie2_summary.xls"
