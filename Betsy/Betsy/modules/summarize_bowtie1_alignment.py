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
        from genomicode import alignlib
        from genomicode import config
        from genomicode import parallel

        log_filenames = _find_output_logs(in_data.identifier)
        assert log_filenames

        results = {}  # dict of sample -> dictionary of output
        for filename in log_filenames:
            # <path>/<sample>.log
            path, file_ = os.path.split(filename)
            f, e = os.path.splitext(file_)
            assert e == ".log"
            sample = f
            results[sample] = alignlib.parse_bowtie1_output(filename)

        # Make table where the rows are the samples and the columns
        # are the statistics.
        all_samples = sorted(results)
        table = []
        header = "Sample", "Aligned Reads", "Total Reads", "Perc Aligned"
        table.append(header)
        for sample in all_samples:
            stats = results[sample]
            total_reads = stats["reads_processed"]
            aligned_reads = stats["aligned_reads"]
            perc_aligned = float(aligned_reads) / total_reads * 100

            x1 = parselib.pretty_int(aligned_reads)
            x2 = parselib.pretty_int(total_reads)
            x3 = "%.2f%%" % perc_aligned
            x = sample, x1, x2, x3
            table.append(x)
        
        # Write out the table as text file.
        TXT_FILE = "summary.txt"
        handle = open(TXT_FILE, 'w')
        for x in table:
            print >>handle, "\t".join(x)
        handle.close()

        txt2xls = filelib.which_assert(config.txt2xls)
        os.system("%s -b %s > %s" % (
            parallel.quote(txt2xls), TXT_FILE, outfile))
            
        
    def name_outfile(self, antecedents, user_options):
        return "bowtie1_summary.xls"


def _find_output_logs(path):
    from genomicode import filelib
    
    # Should be a folder of results from bowtie1.  STDOUT should
    # be sent to ".log" files.
    x = filelib.list_files_in_path(path)
    x = [x for x in x if x.endswith(".log")]
    return x
