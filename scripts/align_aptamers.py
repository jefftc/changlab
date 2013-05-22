#!/usr/bin/env python


DEF_MISMATCH = 0.1
DEF_INSERT = DEF_MISMATCH
DEF_DELBASE = DEF_MISMATCH**2
DEF_DELCHUNK = DEF_MISMATCH**2

OUT_TABLE = "TABLE"
OUT_MARKOV = "MARKOV_MODEL"
OUT_ALIGNMENT = "ALIGNMENT"

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


def _write_table(library, alignment, title, score, is_revcomp, outhandle):
    from genomicode import aptamerlib
    #for x in alignment:
    #    print "ALIGN", repr(x)
    #print "SCORE", score

    # Get the indexes of the regions that are random.
    I_random = [i for i in range(len(library)) if library[i].is_random]
    I_barcode = [i for i in range(len(library)) if library[i].is_barcode]
    assert len(I_random) >= 1, "No random region."
    assert len(I_barcode) <= 1, "Too many barcodes."

    # Count the matches and mismatches.
    counts = {}
    for x in alignment:
        node, log_score, match_type, base_in_lib, base_in_seq = x
        assert node.node_type in [
            aptamerlib.INSERT, aptamerlib.MAIN, aptamerlib.INSERTEND]
        if node.node_type == aptamerlib.INSERTEND:  # not in random region
            continue
        if not library[node.i_seqset].is_random: # only count random region
            continue
        counts[match_type] = counts.get(match_type, 0) + 1
    num_matches = counts.get("MATCH", 0)
    num_mismatches = counts.get("MISMATCH", 0)
    num_deletions = counts.get("DELETE", 0)
    num_insertions = counts.get("INSERT", 0)

    # Pull out the random region.
    ideal_seq, real_seq = aptamerlib.pretty_sequences(alignment)
    ideal_seqs = ideal_seq.split("*")
    actual_seqs = real_seq.split("*")
    assert len(ideal_seqs) >= len(library)
    assert len(actual_seqs) >= len(library)
    x1 = [ideal_seqs[i] for i in I_random]
    #x2 = [actual_seqs[i] for i in I_random]
    ideal_random = "".join(x1)
    #actual_random = "".join(x2)
    x1 = [ideal_seqs[i] for i in I_barcode]
    #x2 = [actual_seqs[i] for i in I_barcode]
    ideal_barcode = "".join(x1)
    #actual_barcode = "".join(x2)
    random_region = ideal_random.replace("-", "")
    barcode = ideal_barcode.replace("-", "")

    # Write the results.
    total_errors = num_mismatches + num_insertions + num_deletions

    # Percent match is hard to calculate.  Should take into
    # account insertions and deletions somehow.
    #perc_match = "%.1f" % (100*float(num_matches) / len(random_region))
    x = title, int(is_revcomp), score, barcode, random_region, \
        len(random_region), total_errors, \
        num_matches, num_mismatches, num_insertions, num_deletions, \
        ideal_seq, real_seq
    print >>outhandle, "\t".join(map(str, x))
    

def _write_alignment(alignment, title, outhandle):
    from genomicode import aptamerlib
    
    ideal_seq, real_seq = aptamerlib.pretty_sequences(alignment)
    print >>outhandle, title
    print >>outhandle, ideal_seq
    print >>outhandle, real_seq
    print >>outhandle


def _write_markov_model(library, alignment, title, sequence, score, is_revcomp,
                        outhandle):
    from genomicode import aptamerlib
    
    print >>outhandle, "%s %s %s" % ("="*10, title, "="*10)
    print >>outhandle, "%s %.2f %d" % (sequence, score, int(is_revcomp))
    for i, x in enumerate(alignment):
        node, log_score, match_type, base_in_lib, base_in_seq = x
        node_str = str(node)
        if node.i_seqset is not None:
            seqset_name = library[node.i_seqset].name
        else:
            assert node.node_type == aptamerlib.INSERTEND
            node_str = "END"
            seqset_name = ""
        x = node_str, seqset_name, match_type, \
            base_in_lib, base_in_seq, "%.2f" % log_score
        print >>outhandle, "\t".join(map(str, x))
    print >>outhandle


def _align_aptamers_h(
    mm, library, base2emission, title, sequence, output, lock):
    import sys
    import StringIO
    from genomicode import aptamerlib
    
    x = aptamerlib.score_sequence(mm, library, base2emission, sequence)
    score, is_revcomp, alignment = x

    outhandle = StringIO.StringIO()
    if output == OUT_TABLE:
        _write_table(library, alignment, title, score, is_revcomp, outhandle)
    elif output == OUT_ALIGNMENT:
        _write_alignment(alignment, title, outhandle)
    elif output == OUT_MARKOV:
        _write_markov_model(
            library, alignment, title, sequence, score, is_revcomp, outhandle)
    else:
        raise AssertionError, "Unknown output format: %s" % output
        
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
        "--output_format", choices=[OUT_TABLE, OUT_MARKOV, OUT_ALIGNMENT],
        default=OUT_TABLE, 
        help="What kind of output to generate.")
    parser.add_argument(
        "--max_alignments", type=int,
        help="Maximum number of sequences to align.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--min_seqlen", type=int, default=None,
        help="Discard sequences less than this minimum length.")
    
    group = parser.add_argument_group(title="Alphabet")
    group.add_argument(
        "--lib2seq", default=[], action="append",
        help="Map the bases in the library to bases in the sequencing.  "
        "If the library contains non-standard bases (i.e. not ACGT), "
        "please indicate how they are read in the sequencing in the format: "
        "<lib_base>:<seq_base>.  E.g. 5:T means that a 5 in the "
        "library is read as a T in the sequence.  If there are "
        "multiple mappings, use this option multiple times.")

    group = parser.add_argument_group(title="Probabilities")
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
    
    args = parser.parse_args()

    assert os.path.exists(args.library_file), \
           "File not found: %s" % args.library_file
    assert os.path.exists(args.sequence_file), \
           "File not found: %s" % args.sequence_file
    assert args.num_procs >= 1 and args.num_procs < 256
    assert args.min_seqlen is None or (
        args.min_seqlen >= 0 and args.min_seqlen < 100)

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

    if args.output_format == OUT_TABLE:
        header = [
            "Title", "Is Revcomp", "Score", "Barcode", "Random Region",
            "Length", "Total Errors", 
            "Num Matches", "Num Mismatches", "Num Insertions", "Num Deletions",
            "Ideal Sequence", "Observed Sequence"]
        print "\t".join(header)
        sys.stdout.flush()   # needed for multiprocessing
    results = []
    for i, x in enumerate(aptamerlib.parse_fastq(args.sequence_file)):
        title, sequence, quality = x

        if args.max_alignments is not None and i >= args.max_alignments:
            break
        if args.min_seqlen is not None and len(sequence) < args.min_seqlen:
            continue
        #if title.find("00013:00070") < 0:
        #    continue

        params = (
            mm, library, base2emission, title, sequence, args.output_format,
            lock)
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
