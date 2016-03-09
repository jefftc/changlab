"""Miscellaneous next generation sequencing tools.

Functions:
get_bedtools_version
make_bedtools_coverage_command
make_bedtools_genomecov_command
parse_bedtools_genomecov_results

"""

def get_bedtools_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel
    
    bedtools = filelib.which_assert(config.bedtools)
    x = parallel.sshell("%s --version" % bedtools, ignore_nonzero_exit=True)
    x = x.strip()
    # bedtools v2.23.0
    # Version: 1.2 (using htslib 1.2.1)
    m = re.search(r"v([\w\. ]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_bedtools_coverage_command(
    bam_filename, features_bed, cov_filename):
    import os
    
    import config
    import filelib
    import parallel

    # XXX Generates a histogram of the counts for each read depth.
    # bedtools coverage [OPTIONS] -abam <align.bam> -b <features.bed>
    bedtools = filelib.which_assert(config.bedtools)
    assert os.path.exists(bam_filename)
    assert os.path.exists(features_bed)

    sq = parallel.quote
    x = [
        sq(bedtools),
        "coverage",
        "-abam", sq(bam_filename),
        "-b", sq(features_bed),
        ">&", sq(cov_filename),
        ]
    return " ".join(x)
    

def make_bedtools_genomecov_command(
    bam_filename, reference_file, cov_filename):
    import os
    import config
    import filelib
    import parallel

    # Generates a histogram of the counts for each read depth.
    # bedtools genomecov [OPTIONS] -ibam <align.bam> -g <ref.fa>
    bedtools = filelib.which_assert(config.bedtools)
    assert os.path.exists(bam_filename)
    assert os.path.exists(reference_file)

    sq = parallel.quote
    x = [
        sq(bedtools),
        "genomecov",
        "-ibam", sq(bam_filename),
        "-g", sq(reference_file),
        ">&", sq(cov_filename),
        ]
    return " ".join(x)


def parse_bedtools_genomecov_results(filename):
    # Return list of (chromosome or "genome", coverage, bases with
    # coverage, size of chromosome, perc bases with coverage).
    # genome  0       18413   4392353 0.00419206
    # genome  1       17191   4392353 0.00391385
    # genome  2       19904   4392353 0.00453151
    # genome  3       27298   4392353 0.00621489
    # ...
    import filelib

    results = []
    for d in filelib.read_row(
        filename, "chrom:s depth:d num_bases:d chr_size:d perc_bases:f"):
        x = d.chrom, d.depth, d.num_bases, d.chr_size, d.perc_bases
        results.append(x)
    return results
