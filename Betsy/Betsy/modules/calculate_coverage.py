from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import ngslib
        from Betsy import module_utils as mlib

        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert filelib.exists_nz(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        features_bed = mlib.get_user_option(
            user_options, "features_bed", check_file=True)
        if features_bed:
            raise NotImplementedError, "features_bed not implemented"

        min_coverage = user_options.get("ignore_coverage_below")
        if min_coverage == "":
            min_coverage = None
        if min_coverage is not None:
            min_coverage = int(min_coverage)
            assert min_coverage >= 0

        metadata["tool"] = "bedtools %s" % ngslib.get_bedtools_version()
        metadata["num_cores"] = num_cores

        # Set up the filenames.
        # list of (sample, bam_filename, cov_filename,
        #   histo_datafile, histo_plotfile, histo_prismfile))
        opj = os.path.join
        jobs = []
        for bam_filename in bam_filenames:
            # <path>/<sample>.bam
            in_path, sample, ext = mlib.splitpath(bam_filename)
            assert ext == ".bam"
            cov_filename = opj(out_path, "%s.coverage.txt" % sample)
            histo_datafile = opj(out_path, "%s.histo.txt" % sample)
            histo_plotfile = opj(out_path, "%s.histo.png" % sample)
            histo_prismfile = opj(out_path, "%s.prism.txt" % sample)
            x = sample, bam_filename, cov_filename, \
                histo_datafile, histo_plotfile, histo_prismfile
            jobs.append(x)

        # Use genomecov to count read depth.
        x = _run_genomecov(jobs, ref_node.identifier, num_cores)
        metadata["commands"] = x

        # Read in the read depth.
        results = {}  # sample -> bedtools results
        for x in jobs:
            sample, bam_filename, cov_filename, \
                    histo_datafile, histo_plotfile, histo_prismfile = x
            assert sample not in results
            results[sample] = ngslib.parse_bedtools_genomecov_results(
                cov_filename)
            
        # Summarize the average read depth.
        _summarize_average_read_depth(results, min_coverage, out_path)

        # Make histograms of the distribution of the read depth for
        # each sample.
        for x in jobs:
            sample, bam_filename, cov_filename, \
                histo_datafile, histo_plotfile, histo_prismfile = x
            _make_histo_file(cov_filename, histo_datafile)
        
        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "coverage"


def _run_genomecov(jobs, reference_file, num_cores):
    from genomicode import parallel
    from genomicode import filelib
    from genomicode import ngslib

    # Set up the commands to run.
    commands = []
    for x in jobs:
        sample, bam_filename, cov_filename, \
                histo_datafile, histo_plotfile, histo_prismfile = x
        x = ngslib.make_bedtools_genomecov_command(
            bam_filename, reference_file, cov_filename)
        commands.append(x)
    parallel.pshell(commands, max_procs=num_cores)

    # Make sure the analysis completed successfully.
    x = [x[2] for x in jobs]
    filelib.assert_exists_nz_many(x)

    return commands


def _summarize_average_read_depth(results, min_coverage, out_path):
    import os
    import math
    from genomicode import parallel
    from genomicode import filelib
    
    # Now calculate the average read depth for each sample.
    sample2depth = {}  # sample -> mean, sd of depth
    for sample in results:
        x = results[sample]
        x = [x for x in x if x[0] == "genome"]
        data = x
        assert data
        # Calculate the mean depth.
        sum_ = 0
        total_bases = 0
        for x in data:
            x, depth, num_bases, x, x = x
            if min_coverage is not None and depth < min_coverage:
                continue
            total_bases += num_bases
            sum_ = sum_ + depth * num_bases
        depth_mean = float(sum_) / total_bases
        # Calculate the SD depth.
        sumsq = 0
        for x in data:
            x, depth, num_bases, x, x = x
            if min_coverage is not None and depth < min_coverage:
                continue
            sumsq += ((depth-depth_mean)*num_bases)**2
        depth_sd = math.sqrt(sumsq / (total_bases-1))
        sample2depth[sample] = depth_mean, depth_sd

    # Make a table of the results.
    mean = "Mean"
    if min_coverage is not None:
        mean = "Mean (>= %d reads)" % min_coverage
    table = []
    header = "Sample", mean, "SD"
    table.append(header)
    for sample in sorted(sample2depth):
        mean, sd = sample2depth[sample]
        mean_s = "%.1f" % mean
        sd_s = "%.2f" % sd
        x = sample, mean_s, sd_s
        table.append(x)

    # Write out the table as text file.
    TXT_FILE = os.path.join(out_path, "summary.txt")
    XLS_FILE = os.path.join(out_path, "summary.xls")
    handle = open(TXT_FILE, 'w')
    for x in table:
        print >>handle, "\t".join(x)
    handle.close()

    sq = parallel.quote
    os.system("txt2xls -b %s > %s" % (sq(TXT_FILE), sq(XLS_FILE)))
    filelib.assert_exists_nz(XLS_FILE)
    os.unlink(TXT_FILE)


def _make_histo_file(cov_filename, outfilename):
    from genomicode import ngslib
    
    results = ngslib.parse_bedtools_genomecov_results(
        cov_filename)

    # [0, 1), [1, 5), [5, 10), etc.
    BINS = [0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 5000, 10000]
    counts = [0] * len(BINS)

    for x in results:
        chrom, coverage, num_bases, chrom_size, perc_bases = x
        if chrom != "genome":
            continue
        bin_i = None
        for i in range(len(BINS)-1):
            b1, b2 = BINS[i], BINS[i+1]
            if coverage >= b1 and coverage < b2:
                bin_i = i
                break
        if bin_i is None:
            assert coverage >= BINS[-1]
            bin_i = len(BINS)-1
        counts[bin_i] += num_bases

    total = sum(counts)
    perc = [float(x)/total for x in counts]
    
    total_noz = sum(counts[1:])
    perc_noz = [float(x)/total_noz for x in counts]
    perc_noz[0] = 0

    total_10 = sum(counts[3:])
    perc_10 = [float(x)/total_10 for x in counts]
    perc_10[0] = perc_10[1] = perc_10[2] = 0

    handle = open(outfilename, 'w')
    header = "Range", "Count", "Fraction", "Fraction (>= 1)", \
             "Fraction (>= 10)"
    print >>handle, "\t".join(header)
    for i in range(len(BINS)):
        if i < len(BINS)-1:
            b1, b2 = BINS[i], BINS[i+1]
            assert b2 > b1
            if b2 == b1+1:
                r = b1
            else:
                r = "%d-%d" % (b1, b2-1)
        else:
            r = ">= %d" % BINS[i]
        x = r, counts[i], perc[i], perc_noz[i], perc_10[i]
        assert len(x) == len(header)
        print >>handle, "\t".join(map(str, x))
    handle.close()
    

