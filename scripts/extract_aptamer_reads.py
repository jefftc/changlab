#!/usr/bin/env python

def main():
    import os
    import sys
    import time
    import argparse
    #import multiprocessing

    from genomicode import aptamerlib
    from genomicode import parselib

    parser = argparse.ArgumentParser(
        description="Pull out the reads that originate from an aptamer "
        "library.")
    parser.add_argument("library_file")
    parser.add_argument(
        "sequence_file", help="FASTQ-formatted sequence file.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--min_seqlen", type=int, default=None,
        help="Discard sequences less than this minimum length.")

    parser.add_argument(
        "--match_file", help="File for reads that match this library.")
    parser.add_argument(
        "--leftover_file", help="File for leftover reads that don't match.")
    parser.add_argument("--clobber", default=False, action="store_true")

    args = parser.parse_args()

    # Check the inputs.
    assert os.path.exists(args.library_file), \
           "File not found: %s" % args.library_file
    assert os.path.exists(args.sequence_file), \
           "File not found: %s" % args.sequence_file
    assert args.num_procs >= 1 and args.num_procs < 256
    assert args.min_seqlen is None or (
        args.min_seqlen >= 0 and args.min_seqlen < 100)
    
    assert args.match_file, "Please specify a match_file."
    if not args.clobber and (
        args.match_file and os.path.exists(args.match_file)):
        raise AssertionError, ("match_file %s exists.  "
              "Please use --clobber to overwrite." % args.match_file)
    if not args.clobber and (
        args.leftover_file and os.path.exists(args.leftover_file)):
        raise AssertionError, ("leftover_file %s exists.  "
              "Please use --clobber to overwrite." % args.leftover_file)

    match_handle = open(args.match_file, 'w')
    leftover_handle = None
    if args.leftover_file:
        leftover_handle = open(args.leftover_file, 'w')

    library = aptamerlib.read_library(args.library_file)

    #manager = multiprocessing.Manager()
    #lock = manager.Lock()
    #pool = multiprocessing.Pool(args.num_procs)

    TIME_FORMAT = "%m/%d/%Y %H:%M:%S"
    last_time = None
    for i, x in enumerate(aptamerlib.parse_fastq(args.sequence_file)):
        title, sequence, quality = x

        if args.min_seqlen is not None and len(sequence) < args.min_seqlen:
            continue

        t = time.time()
        if last_time is None or t > last_time + 5:
            last_time = t
            now = time.strftime(TIME_FORMAT, time.localtime(t))
            print "%s\t%s" % (
                now, "Extracting read %s." % parselib.pretty_int(i+1))
            sys.stdout.flush()

        orientation = aptamerlib.guess_sequence_orientation(sequence, library)
        assert orientation in [-1, 0, 1]

        if orientation in [-1, 1]:
            # Matches the library.
            print >>match_handle, title
            print >>match_handle, sequence
            print >>match_handle, "+"
            print >>match_handle, quality
        else:
            print >>leftover_handle, title
            print >>leftover_handle, sequence
            print >>leftover_handle, "+"
            print >>leftover_handle, quality
            
    match_handle.close()
    leftover_handle.close()
            

if __name__ == '__main__':
    main()
