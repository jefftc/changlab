#!/usr/bin/env python


DEF_MISMATCH = 0.1
DEF_INSERT = DEF_MISMATCH
DEF_DELBASE = DEF_MISMATCH**2
DEF_DELCHUNK = DEF_MISMATCH**2



def parse_fastq(filename):
    # Iterator that yields tuples (title, sequence, quality).
    from genomicode import filelib

    # Format of FASTQ files:
    # @4GEOU:00042:00049                          Title
    # ACTGCTAATTCACACTGGATTAGTTGGGCTACTTCATCGT    Sequence
    # +                                           Always "+"
    # =<>>A77.7.54>4444.46-444,44*3333:9:44443    Quality
    
    handle = filelib.openfh(filename)
    while True:
        x = [handle.readline() for i in range(4)]
        lines = [x.strip() for x in x]
        if not lines[0]:
            break
        title, sequence, x, quality = lines
        assert x == "+"
        assert len(sequence) == len(quality)
        assert quality

        yield title, sequence, quality

def _parse_base2emission(base2emission_list):
    # base2emission is a list of "<lib_base>:<seq_base>".  Return a
    # dictionary of <lib_base> -> <seq_base>.
    base2emission = {}
    for x in base2emission_list:
        assert ":" in x, "improper format: %s" % x
        y = x.split(":")
        assert len(y) == 2, "improper format: %s" % x
        lib_base, seq_base = y
        assert lib_base not in base2emission, "multiple: %s" % lib_base
        base2emission[lib_base] = seq_base
    return base2emission
    

def _align_aptamers_h(mm, library, base2emission, title, sequence, lock):
    import sys
    import StringIO
    from genomicode import aptamerlib
    
    x = aptamerlib.score_sequence(mm, library, base2emission, sequence)
    score, is_revcomp, alignment = x

    # Get the indexes of the regions that are random.
    I_random = [i for (i, (name, is_random, seqset)) in enumerate(library)
                if is_random]

    outhandle = StringIO.StringIO()
    if 1:
        ideal_seq, real_seq = aptamerlib.pretty_sequences(alignment)

        ideal_seqs = ideal_seq.split("*")
        actual_seqs = real_seq.split("*")
        assert len(ideal_seqs) >= len(library)
        assert len(actual_seqs) >= len(library)

        # Pull out the random region.
        x1 = [ideal_seqs[i] for i in I_random]
        x2 = [actual_seqs[i] for i in I_random]
        ideal_random = "".join(x1)
        actual_random = "".join(x2)

        # Count the matches and mismatches.
        assert len(ideal_random) == len(actual_random)
        num_matches = num_mismatches = num_insertions = num_deletions = 0
        for (ideal, actual) in zip(ideal_random, actual_random):
            assert not (ideal == "-" and actual == "-")
            if base2emission.get(ideal, ideal) == actual:
                num_matches += 1
            elif actual == "-":
                num_deletions += 1
            elif ideal == "-":
                num_insertions += 1
            else:
                num_mismatches += 1
        
        # Write the results.
        random_region = ideal_random.replace("-", "")
        total_errors = num_mismatches + num_insertions + num_deletions
        perc_match = "%.1f" % (100*float(num_matches) / len(random_region))
        x = title, int(is_revcomp), score, random_region, \
            len(random_region), perc_match, total_errors, \
            num_matches, num_mismatches, num_insertions, num_deletions, \
            ideal_seq, real_seq
        print >>outhandle, "\t".join(map(str, x))
    else:
        # Alternate output format.
        print >>outhandle, "%s %s %s" % ("="*10, title, "="*10)
        print >>outhandle, sequence
        print >>outhandle, score, is_revcomp
        for i, x in enumerate(alignment):
            node, match_type, base_in_lib, base_in_seq = x
            node_str = aptamerlib.format_node(node)
            if len(node) >= 4:
                x, i_seqset, i_sequence, i_base = node
                seqset_name = library[i_seqset][0]
            else:
                assert node[0] == aptamerlib.INSERTEND
                node_str = "END"
                seqset_name = ""
            x = node_str, seqset_name, match_type, \
                base_in_lib, base_in_seq
            print >>outhandle, "\t".join(map(str, x))
        print >>outhandle
        
    lock.acquire()
    sys.stdout.write(outhandle.getvalue())
    sys.stdout.flush()
    lock.release()


def main():
    import os
    import sys
    import argparse
    import multiprocessing

    from genomicode import aptamerlib

    parser = argparse.ArgumentParser(
        description="Align sequences to an aptamer library.")
    parser.add_argument("library_file")
    parser.add_argument(
        "sequence_file", help="FASTQ-formatted sequence file.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--max_alignments", type=int,
        help="Maximum number of sequences to align.")

    group = parser.add_argument_group(title="Model Specification")
    group.add_argument(
        "--mismatch", type=float, default=DEF_MISMATCH,
        help="Probability of a mismatch (default %0.2f)." % DEF_MISMATCH)
    group.add_argument(
        "--insert", type=float, default=DEF_INSERT,
        help="Probability of an insertion (default %0.2f)." % DEF_INSERT)
    group.add_argument(
        "--delete_base", type=float, default=DEF_DELBASE,
        help="Probability of a deletion (default %0.2f)." % DEF_DELBASE)
    group.add_argument(
        "--delete_chunk", type=float, default=DEF_DELCHUNK,
        help="Probability of a deletion (default %0.2f)." % DEF_DELCHUNK)
    group.add_argument(
        "--lib2seq", default=[], action="append",
        help="Map the bases in the library to bases in the sequencing.  "
        "If the library contains non-standard bases (i.e. not ACGT), "
        "please indicate how they are read in the sequencing in the format: "
        "<lib_base>:<seq_base>.  E.g. 5:T means that a 5 in the "
        "library is read as a T in the sequence.  If there are "
        "multiple mappings, use this option multiple times.")
    
    args = parser.parse_args()

    assert os.path.exists(args.library_file), \
           "File not found: %s" % args.library_file
    assert os.path.exists(args.sequence_file), \
           "File not found: %s" % args.sequence_file
    assert args.num_procs >= 1 and args.num_procs < 256

    assert args.mismatch >= 0 and args.mismatch < 0.5, \
           "mismatch (%s) should be between 0 and 0.5" % args.mismatch
    assert args.insert >= 0 and args.insert < 0.5, \
           "insert (%s) should be between 0 and 0.5" % args.insert
    assert args.delete_base >= 0 and args.delete_base < 0.5, \
           "delete_base (%s) should be between 0 and 0.5" % args.delete_base
    assert args.delete_chunk >= 0 and args.delete_chunk < 0.5, \
           "delete_chunk (%s) should be between 0 and 0.5" % args.delete_chunk
    assert args.insert + args.delete_base + args.delete_chunk <= 0.5, \
           "Too much probability allocated for indels."

    base2emission = _parse_base2emission(args.lib2seq)

    library = aptamerlib.read_library(args.library_file)
    mm = aptamerlib.make_markov_model(
        library, base2emission, p_mismatch=args.mismatch, p_insert=args.insert,
        p_delbase=args.delete_base, p_delchunk=args.delete_chunk)


    manager = multiprocessing.Manager()
    lock = manager.Lock()
    pool = multiprocessing.Pool(args.num_procs)

    header = [
        "Title", "Is Revcomp", "Score", "Random Region",
        "Length", "Percent Aligned", "Total Errors", 
        "Num Matches", "Num Mismatches", "Num Insertions", "Num Deletions",
        "Ideal Sequence", "Observed Sequence"]
    print "\t".join(header)
    sys.stdout.flush()
    results = []
    for i, x in enumerate(parse_fastq(args.sequence_file)):
        title, sequence, quality = x

        if args.max_alignments is not None and i >= args.max_alignments:
            break

        params = mm, library, base2emission, title, sequence, lock
        keywds = {}
        if args.num_procs == 1:
            _align_aptamers_h(*params, **keywds)
        else:
            x = pool.apply_async(_align_aptamers_h, params, keywds)
            results.append(x)
    pool.close()
    pool.join()
    for x in results:
        x.get()


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
