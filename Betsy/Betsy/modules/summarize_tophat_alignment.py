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

        align_node = in_data
        x = filelib.list_files_in_path(
            align_node.identifier, endswith="align_summary.txt")
        align_filenames = x
        assert align_filenames, "Missing align_summary.txt"

        results = {}  # dict of sample -> dictionary of output
        for filename in align_filenames:
            # Names must in the format:
            # <path>/<sample>.tophat/alignment_summary.txt
            # full_path   <path>/<sample>.tophat
            # path        <path>
            # tophat_dir  <sample>.tophat
            # file_       accepted_hits.bam
            # sample      <sample>
            
            full_path, file_ = os.path.split(filename)
            path, tophat_dir = os.path.split(full_path)
            assert file_ == "align_summary.txt"
            assert tophat_dir.endswith(".tophat")
            sample = tophat_dir[:-7]

            x = alignlib.parse_tophat_align_summary(filename)
            results[sample] = x

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
        return "tophat_summary.xls"
