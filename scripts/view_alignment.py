#!/usr/bin/env python

def _parse_position(pos_str):
    # Return chrom, pos (as integer)
    # chr20:45,927,663
    # chr20:45927663
    assert ":" in pos_str, "No colon: %s" % pos_str
    chrom, pos = pos_str.split(":", 2)
    pos = pos.replace(",", "")
    assert "-" not in pos, "Invalid position: %s" % pos
    pos = int(pos)
    assert pos >= 0
    return chrom, pos


def main():
    import os
    import argparse
    import itertools
    
    from genomicode import filelib
    from genomicode import config
    from genomicode import parallel

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("reference_genome", help="fastq file")

    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--dry_run", action="store_true",
        help="Just display the commands, and don't generate the alignment.")
    parser.add_argument(
        "--window", default=80, type=int,
        help="Number of bases in alignment.  Default: 80")

    group = parser.add_argument_group(title="Input")
    group.add_argument(
        "--bam_file", help="Indexed BAM file.")
    group.add_argument(
        "--bam_path", help="Path to BAM files.")
    group.add_argument(
        "--position", action="append", default=[],
        help="Specify a position to view, "
        "e.g. chr20:45,927,663 or chr20:45927663.  1-based coordinates")
    group.add_argument(
        "--position_file", 
        help="Tab-delimited text file with two columns.  "
        "Column 1 is chromosome, column 2 is position.")
    
    group = parser.add_argument_group(title="Output")
    group.add_argument(
        "--prefix", help="Pre-pend a prefix to each outfile.")
    group.add_argument(
        "--outpath", help="If multiple alignments are generated, this option "
        "directs where to save the output files.")


    # Parse the input arguments.
    args = parser.parse_args()
    filelib.assert_exists_nz(args.reference_genome)
    assert args.bam_file or args.bam_path, \
           "Either --bam_file or --bam_path must be provided."
    assert not (args.bam_file and args.bam_path), \
           "Cannot specify both --bam_file or --bam_path."
    if args.bam_file:
        filelib.assert_exists_nz(args.bam_file)
    if args.bam_path:
        assert os.path.exists(args.bam_path)
    if args.position_file:
        filelib.assert_exists_nz(args.position_file)
    if args.outpath and not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    if args.num_procs < 1 or args.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")
    assert args.window >= 1 and args.window < 500

    bam_filenames = []
    if args.bam_file:
        bam_filenames.append(args.bam_file)
    else:
        x = os.listdir(args.bam_path)
        x = [x for x in x if x.endswith(".bam")]
        x = [os.path.join(args.bam_path, x) for x in x]
        bam_filenames = x
    assert bam_filenames, "No bam files found."

    positions = []  # list of (chrom, pos)
    for x in args.position:
        chrom, pos = _parse_position(x)
        positions.append((chrom, pos))
    if os.path.exists(args.position_file):
        for cols in filelib.read_cols(args.position_file):
            assert len(cols) == 2
            chrom, pos = cols
            pos = int(pos)
            assert pos >= 1
            positions.append((chrom, pos))
    assert positions, "No positions specified."

    # Make the commands.
    assert hasattr(config, "samtools")
    filelib.assert_exists(config.samtools)

    commands = []
    for x in itertools.product(bam_filenames, positions):
        bam_filename, (chrom, pos) = x
        
        p, f = os.path.split(bam_filename)
        sample, e = os.path.splitext(f)

        left = max(pos-args.window/2, 1)
        pos_str = "%s:%s" % (chrom, left)

        x = "%2s.%9s.%s.html" % (chrom, pos, sample)
        if args.prefix:
            x = "%s.%s" % (args.prefix, x)
        if args.outpath:
            x = os.path.join(args.outpath, x)
        out_filename = x

        # samtools tview -d t -p 7:100550778 bam01/196B-lung.bam $FA
        sq = parallel.quote
        x = [
            sq(config.samtools), "tview",
            "-d", "h",
            "-p", pos_str,
            sq(bam_filename),
            sq(args.reference_genome),
            ]
        x = " ".join(x)
        x = "%s >& %s" % (x, sq(out_filename))
        commands.append(x)

    if args.dry_run:
        for x in commands:
            print x
        return

    parallel.pshell(commands, max_procs=args.num_procs)
    

if __name__ == '__main__':
    main()
