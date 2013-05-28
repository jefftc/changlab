#!/usr/bin/env python

def _parse_titles(titles):
    # ["title1,title2", "title3"]
    all_titles = []
    for x in titles:
        x = x.split(",")
        all_titles.extend(x)
    return all_titles

def main():
    import os
    import sys
    import time
    import argparse
    #import multiprocessing

    from genomicode import aptamerlib
    from genomicode import parselib

    parser = argparse.ArgumentParser(description="Pull out a subset of reads.")
    parser.add_argument(
        "sequence_file", help="FASTQ-formatted sequence file.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")

    parser.add_argument(
        "--match_file", help="File for reads that match this library.")
    parser.add_argument(
        "--leftover_file", help="File for leftover reads that don't match.")
    parser.add_argument("--clobber", default=False, action="store_true")

    parser.add_argument(
        "--min_seqlen", type=int, default=None,
        help="Discard sequences less than this minimum length.")
    parser.add_argument(
        "--library_file", help="Want reads that match this library.")
    parser.add_argument(
        "--titles", default=[], action="append",
        help="Want reads with these titles.  "
        "Comma-separated titles, parameter can be used multiple times.")

    args = parser.parse_args()

    # Check the inputs.
    assert os.path.exists(args.sequence_file), \
           "File not found: %s" % args.sequence_file
    assert args.num_procs >= 1 and args.num_procs < 256
    assert args.min_seqlen is None or (
        args.min_seqlen >= 0 and args.min_seqlen < 100)
    assert not args.library_file or os.path.exists(args.library_file), \
           "File not found: %s" % args.library_file
    
    assert args.match_file, "Please specify a match_file."
    if not args.clobber and (
        args.match_file and os.path.exists(args.match_file)):
        raise AssertionError, ("match_file %s exists.  "
              "Please use --clobber to overwrite." % args.match_file)
    if not args.clobber and (
        args.leftover_file and os.path.exists(args.leftover_file)):
        raise AssertionError, ("leftover_file %s exists.  "
              "Please use --clobber to overwrite." % args.leftover_file)

    titles = _parse_titles(args.titles)

    match_handle = open(args.match_file, 'w')
    leftover_handle = None
    if args.leftover_file:
        leftover_handle = open(args.leftover_file, 'w')

    library = None
    if args.library_file:
        library = aptamerlib.read_library(args.library_file)

    #manager = multiprocessing.Manager()
    #lock = manager.Lock()
    #pool = multiprocessing.Pool(args.num_procs)

    TIME_FORMAT = "%m/%d/%Y %H:%M:%S"
    last_time = None
    for i, x in enumerate(aptamerlib.parse_fastq(args.sequence_file)):
        title, sequence, quality = x

        t = time.time()
        if last_time is None or t > last_time + 5:
            last_time = t
            now = time.strftime(TIME_FORMAT, time.localtime(t))
            print "%s\t%s" % (
                now, "Extracting read %s." % parselib.pretty_int(i+1))
            sys.stdout.flush()

        if args.min_seqlen is not None and len(sequence) < args.min_seqlen:
            continue

        is_match = False
        
        # Keep if either the title matches or the library matches.
        if titles:
            if title in titles:
                is_match = True
        if library:
            orientation = aptamerlib.guess_sequence_orientation(
                sequence, library)
            assert orientation in [-1, 0, 1]
            if orientation in [-1, 1]:
                is_match = True

        if is_match:
            # Matches the library.
            print >>match_handle, title
            print >>match_handle, sequence
            print >>match_handle, "+"
            print >>match_handle, quality
        elif leftover_handle:
            print >>leftover_handle, title
            print >>leftover_handle, sequence
            print >>leftover_handle, "+"
            print >>leftover_handle, quality
            
    match_handle.close()
    if leftover_handle:
        leftover_handle.close()
            

if __name__ == '__main__':
    main()
