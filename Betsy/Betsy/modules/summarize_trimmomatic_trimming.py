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
        from genomicode import config
        from genomicode import parallel

        log_filenames = find_trimmomatic_logs(in_data.identifier)
        assert log_filenames

        results = {}  # dict of sample -> dictionary of output
        for filename in log_filenames:
            # <path>/<sample>.log
            path, file_ = os.path.split(filename)
            f, e = os.path.splitext(file_)
            assert e == ".log"
            sample = f
            results[sample] = parse_trimmomatic_output(filename)

        # Make table where the rows are the samples and the columns
        # are the statistics.
        all_samples = sorted(results)
        table = []
        header = "Sample", "Total Reads", "Dropped Reads", "Perc Dropped"
        table.append(header)
        for sample in all_samples:
            stats = results[sample]
            reads = stats["reads_processed"]
            dropped = stats["dropped_reads"]
            perc = float(dropped) / reads * 100

            x1 = parselib.pretty_int(reads)
            x2 = parselib.pretty_int(dropped)
            x3 = "%.2f%%" % perc
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
        return "trimmomatic_summary.xls"


def find_trimmomatic_logs(path):
    from genomicode import filelib
    
    # Should be a folder of results from trimmomatic.  STDOUT should
    # be sent to ".log" files.
    x = filelib.list_files_in_path(path)
    x = [x for x in x if x.endswith(".log")]
    return x

def parse_trimmomatic_output(filename):
    # Return a dictionary with keys:
    # reads_processed
    # dropped_reads
    #
    # For paired ends, this refers to pairs of reads.  dropped_reads
    # indicates how many pairs were dropped, because on or both were
    # dropped.
    from genomicode import filelib

    # Single end reads:
    # Input Reads: 1764254 Surviving: 1764160 (99.99%) Dropped: 94 (0.01%)
    #
    # Paired end reads:
    # Input Read Pairs: 60032406 Both Surviving: 59093198 (98.44%)
    #   Forward Only Surviving: 891164 (1.48%) Reverse Only Surviving:
    #   20511 (0.03%) Dropped: 27533 (0.05%)

    results = {}
    for line in filelib.openfh(filename):
        if not line.startswith("Input Read"):
            continue
        cols = line.strip().split()
        if line.startswith("Input reads:"):
            # Single end.
            assert len(cols) == 9
            reads = int(cols[2])
            dropped = int(cols[7])
        elif line.startswith("Input Read Pairs:"):
            # Paired end.
            assert len(cols) == 21
            reads = int(cols[3])
            dropped = reads - int(cols[6])
        assert dropped < reads
        results["reads_processed"] = reads
        results["dropped_reads"] = dropped
    assert results, "Parse error: %s" % filename
    return results

