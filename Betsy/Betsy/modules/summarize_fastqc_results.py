from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os

        # Should be a folder of fastqc results.
        in_path = in_data.identifier
        x = os.listdir(in_path)
        x = [x for x in x if x.endswith("_fastqc")]
        x = [os.path.join(in_path, x, "summary.txt") for x in x]
        x = [x for x in x if os.path.exists(x)]
        filenames = x
        assert filenames, "No FastQC results found."

        results = []
        for filename in filenames:
            x = _read_fastqc_summary(filename)
            results.extend(x)

        # Clean up the filenames (to sample names)
        for i in range(len(results)):
            status, statistic, sample = results[i]
            if sample.lower().endswith(".gz"):
                sample = sample[:-3]
            if sample.lower().endswith(".fq"):
                sample = sample[:-3]
            if sample.lower().endswith(".fastq"):
                sample = sample[:-6]
            results[i] = status, statistic, sample

        # Make table where the rows are the samples and the columns
        # are the statistics.
        all_statistics = []  # statistics, preserve order
        for x in results:
            status, statistic, sample = x
            if statistic not in all_statistics:
                all_statistics.append(statistic)
        # Make list of samples, alphabetical.
        x = [x[-1] for x in results]
        x = sorted({}.fromkeys(x))
        all_samples = x

        sample2stat2status = {}
        for x in results:
            status, statistic, sample = x
            if sample not in sample2stat2status:
                sample2stat2status[sample] = {}
            sample2stat2status[sample][statistic] = status

        table = []
        for sample in all_samples:
            x = [sample2stat2status.get(sample, {}).get(x, "")
                 for x in all_statistics]
            table.append(x)

        # Add the sample name to the first column of the table.
        assert len(all_samples) == len(table)
        for i in range(len(table)):
            table[i] = [all_samples[i]] + table[i]
        # Add the statistic name to the first row.
        x = ["FastQC"] + all_statistics
        table = [x] + table

        # Write out the table as text file.
        TXT_FILE = "fastqc_summary.txt"
        handle = open(TXT_FILE, 'w')
        for x in table:
            print >>handle, "\t".join(x)
        handle.close()

        os.system("txt2xls -b %s > %s" % (TXT_FILE, outfile))
            
        
    def name_outfile(self, antecedents, user_options):
        return "fastqc_summary.xls"

    
def _read_fastqc_summary(filename):
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
