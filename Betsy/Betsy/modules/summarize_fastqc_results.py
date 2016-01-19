from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib
        from genomicode import sortlib

        # Should be a folder of fastqc results.
        fastqc_path = in_data.identifier

        # Find all the FASTQC results.
        x = filelib.list_files_in_path(fastqc_path, endswith="summary.txt")
        x = [os.path.split(x)[0] for x in x]
        paths = x
        assert paths, "No FASTQC files found."

        # Read the results.
        all_results = [read_fastqc_results(x) for x in paths]
        assert all_results

        # Make table where the rows are the samples and the columns
        # are the statistics.
        sample2results = {}
        for x in all_results:
            assert x.sample not in sample2results
            sample2results[x.sample] = x
        all_statistics = all_results[0].statistics_order
        all_samples = sortlib.sort_natural(sample2results)

        table = []
        header = [
            "Sample", "Total Sequences", "Filtered Sequences",
            "Sequence length", "GC"] + all_statistics
        table.append(header)
        for sample in all_samples:
            results = sample2results[sample]
            x1 = [sample]
            x2 = [
                results.total_sequences, results.filtered_sequences,
                results.sequence_length, results.percent_gc]
            x3 = [results.statistics[x] for x in all_statistics]
            x = x1 + x2 + x3
            assert len(x) == len(header)
            table.append(x)

        # Write out the table as text file.
        TXT_FILE = "fastqc_summary.txt"
        handle = open(TXT_FILE, 'w')
        for x in table:
            print >>handle, "\t".join(map(str, x))
        handle.close()

        os.system("txt2xls -b %s > %s" % (TXT_FILE, outfile))
            
        
    def name_outfile(self, antecedents, user_options):
        return "fastqc_summary.xls"


class FastQCResults:
    def __init__(self, sample, total_sequences, filtered_sequences,
                 sequence_length, percent_gc, statistics, statistics_order):
        # statistics is a dictionary of name of statistic -> status
        # statistics_order is the order that the statistics were given
        # in the fastqc output.
        assert sorted(statistics) == sorted(statistics_order)
        self.sample = sample
        self.total_sequences = total_sequences
        self.filtered_sequences = filtered_sequences
        self.sequence_length = sequence_length
        self.percent_gc = percent_gc
        self.statistics = statistics.copy()
        self.statistics_order = statistics_order[:]


def read_fastqc_results(fastqc_path):
    import os
    from genomicode import filelib

    summary_file = os.path.join(fastqc_path, "summary.txt")
    data_file = os.path.join(fastqc_path, "fastqc_data.txt")
    filelib.assert_exists_nz(summary_file)
    filelib.assert_exists_nz(data_file)

    summary = read_fastqc_summary(summary_file)
    data = read_fastqc_data(data_file)

    # Figure out the sample names from the filenames.
    samples = sorted([x[-1] for x in summary])
    assert samples[0] == samples[-1], "%s %s" % (samples[0], samples[-1])
    sample = samples[0]
    if sample.lower().endswith(".gz"):
        sample = sample[:-3]
    if sample.lower().endswith(".fq"):
        sample = sample[:-3]
    if sample.lower().endswith(".fastq"):
        sample = sample[:-6]

    # Make the statistics dictionary.
    statistics = {}
    statistics_order = []
    for x in summary:
        status, statistic, x = x
        assert statistic not in statistics
        statistics[statistic] = status
        statistics_order.append(statistic)

    x = FastQCResults(
        sample, data["total_sequences"], data["filtered_sequences"],
        data["sequence_length"], data["percent_gc"],
        statistics, statistics_order)
    return x
    

def read_fastqc_summary(filename):
    # Return list of (<status>, <statistic>, <filename>)
    import os
    from genomicode import filelib
    
    assert os.path.exists(filename)
    data = []
    for x in filelib.read_cols(filename):
        assert len(x) == 3
        status, statistic, filename = x
        data.append((status, statistic, filename))
    return data


def read_fastqc_data(filename):
    # Return a dictionary of:
    # total_sequences     <int>
    # filtered_sequences  <int>
    # sequence_length     <str>    "205", "15-205"
    # percent_gc          <float>

    data = {}
    for line in open(filename):
        # Line seems to end with:
        # 'Total Sequences\t1056547\t\n'
        # Not enough just to strip \r\n.
        #cols = line.rstrip("\r\n").split("\t")
        cols = line.rstrip().split("\t")

        if line.startswith("Total Sequences"):
            assert len(cols) == 2, repr(line)
            data["total_sequences"] = int(cols[1])
        elif line.startswith("Filtered Sequences"):
            assert len(cols) == 2
            data["filtered_sequences"] = int(cols[1])
        elif line.startswith("Sequence length"):
            assert len(cols) == 2
            data["sequence_length"] = cols[1]
        elif line.startswith("%GC"):
            assert len(cols) == 2
            data["percent_gc"] = float(cols[1])/100
    assert len(data) == 4, "Error parsing: %s" % filename
    return data
