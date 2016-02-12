from Module import AbstractModule


class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import math
        from genomicode import shell
        from genomicode import filelib
        from genomicode import ngslib
        from Betsy import module_utils

        bam_node, reference_node = antecedents
        bam_filenames = filelib.list_files_in_path(
            bam_node.identifier, endswith=".bam")
        reference_file = reference_node.identifier
        filelib.safe_mkdir(out_path)

        assert os.path.exists(reference_file)

        features_bed = module_utils.get_user_option(
            user_options, "features_bed", check_file=True)
        if features_bed:
            raise NotImplementedError, "features_bed not implemented"

        # Set up bedtools genomecov jobs.
        jobs = []   # list of sample, bam_filename, coverage_filename
        for bam_filename in bam_filenames:
            # <path>/<sample>.bam
            in_path, file_ = os.path.split(bam_filename)
            f, e = os.path.splitext(file_)
            assert e == ".bam"
            sample = f
            coverage_filename = os.path.join(
                out_path, "%s.coverage.txt" % sample)
            jobs.append((sample, bam_filename, coverage_filename))

        # Set up the commands to run.
        commands = []
        for x in jobs:
            sample, bam_filename, cov_filename = x
            x = ngslib.make_bedtools_genomecov_command(
                bam_filename, reference_file, cov_filename)
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)
        
        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, bam_filename, cov_filename = x
            assert filelib.exists_nz(cov_filename), \
                   "Missing: %s" % cov_filename

        results = {}  # sample -> bedtools results
        for x in jobs:
            sample, bam_filename, cov_filename = x
            assert sample not in results
            results[sample] = ngslib.parse_bedtools_genomecov_results(
                cov_filename)

        min_coverage = user_options.get("ignore_coverage_below")
        if min_coverage == "":
            min_coverage = None
        if min_coverage is not None:
            min_coverage = int(min_coverage)
            assert min_coverage >= 0
        
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

        sq = shell.quote
        os.system("txt2xls -b %s > %s" % (sq(TXT_FILE), sq(XLS_FILE)))
            
        # TODO: 
        # Make plots showing the histogram of the read depths.

        
    def name_outfile(self, antecedents, user_options):
        return "trimmomatic_summary.xls"
