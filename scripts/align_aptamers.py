#!/usr/bin/env python


DEF_MATCH = 0.9
DEF_MISMATCH = (1.0-DEF_MATCH)

DEF_NOINDEL = 0.9
DEF_INSERT = (1.0-DEF_NOINDEL)
DEF_DELETE = (1.0-DEF_NOINDEL)

# Have different probabilities for INDELs in random region.  Lower
# probability means less likely to have insertions and deletions here.
DEF_INSERT_RANDOM = DEF_INSERT/2.0
DEF_DELETE_RANDOM = DEF_DELETE/2.0   # NOT USED


#OUT_TABLE = "TABLE"
#OUT_MARKOV = "MARKOV_MODEL"
#OUT_ALIGNMENT = "ALIGNMENT"


# _parse_titles
# _parse_base2emission
#
# _write_table
# _write_alignment
# _write_markov_format
#
# _align_aptamers_h_h
# _align_aptamers_h
# _append_to_file_or_handle
#
# write_markov_model


def _parse_titles(titles):
    # ["title1,title2", "title3"]
    all_titles = []
    for x in titles:
        x = x.split(",")
        all_titles.extend(x)
    return all_titles

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
    x2 = [actual_seqs[i] for i in I_random]
    ideal_random = "".join(x1)
    actual_random = "".join(x2)
    x1 = [ideal_seqs[i] for i in I_barcode]
    x2 = [actual_seqs[i] for i in I_barcode]
    ideal_barcode = "".join(x1)
    actual_barcode = "".join(x2)
    #random_region = ideal_random.replace("-", "")
    random_region = actual_random
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


def _write_markov_format(
    library, alignment, title, sequence, score, is_revcomp, outhandle):
    import math
    from genomicode import aptamerlib
    
    print >>outhandle, "%s %s %s" % ("="*10, title, "="*10)
    print >>outhandle, "%s %.2f %d" % (sequence, score, int(is_revcomp))
    prev_score = 0.0
    for i, x in enumerate(alignment):
        node, log_score, match_type, base_in_lib, base_in_seq = x
        node_str = str(node)
        
        score_diff = 0.0
        if log_score < prev_score:
            log_score_diff = log_score - prev_score
            prev_score = log_score
            score_diff = math.exp(log_score_diff)
        score_diff = "%.3g" % score_diff
        
        if node.i_seqset is not None:
            seqset_name = library[node.i_seqset].name
        else:
            assert node.node_type == aptamerlib.INSERTEND
            node_str = "END"
            seqset_name = ""
        x = i+1, node_str, seqset_name, match_type, \
            base_in_lib, base_in_seq, "%.2f" % log_score, score_diff
        print >>outhandle, "\t".join(map(str, x))
    print >>outhandle


def _align_aptamers_h_h(
    mm, library, base2emission, title, sequence, 
    table_handle, alignment_handle, markov_handle):
    from genomicode import aptamerlib
    
    x = aptamerlib.score_sequence(mm, library, base2emission, sequence)
    score, is_revcomp, alignment = x

    assert table_handle
    _write_table(library, alignment, title, score, is_revcomp, table_handle)
    if alignment_handle:
        _write_alignment(alignment, title, alignment_handle)
    if markov_handle:
        _write_markov_format(
            library, alignment, title, sequence, score, is_revcomp,
            markov_handle)


def _append_to_file_or_handle(file_or_handle, s):
    if file_or_handle is None:
        return
    handle = file_or_handle
    if type(file_or_handle) is type(""):
        handle = open(file_or_handle, 'a')
    handle.write(s)
    handle.flush()
    if type(file_or_handle) is type(""):
        handle.close()
    
        

def _align_aptamers_h(
    start, skip, sequence_file, max_alignments, min_seqlen, titles, 
    mm, library, base2emission, alignment_file, markov_file,
    lock):
    import sys
    from StringIO import StringIO
    from genomicode import aptamerlib

    table_handle = sys.stdout   # can't pass handle to multiproc function

    table_cache = StringIO()
    alignment_cache = StringIO()
    markov_cache = StringIO()

    files = [
        (table_handle, table_cache),
        (alignment_file, alignment_cache),
        (markov_file, markov_cache),
        ]
    files = [x for x in files if x[0]]  # no missing handles
    
    for i, x in enumerate(aptamerlib.parse_fastq(sequence_file)):
        title, sequence, quality = x
        if i % skip != start:
            continue
        if max_alignments is not None and i >= max_alignments:
            break
        if min_seqlen is not None and len(sequence) < min_seqlen:
            continue
        if titles and title not in titles:
            continue

        _align_aptamers_h_h(
            mm, library, base2emission, title, sequence, 
            table_cache, alignment_cache, markov_cache)

        # ~266 bytes per alignment.  Markov output is ~17x size of
        # table output.
        if len(table_cache.getvalue()) >= 256*1024:
            lock.acquire()
            for h, c in files:
                _append_to_file_or_handle(h, c.getvalue())
                c.truncate(0)
            lock.release()
    lock.acquire()
    for h, c in files:
        _append_to_file_or_handle(h, c.getvalue())
    lock.release()


def write_markov_model(library, base2emission, p_mismatch, p_insert, p_delete,
                       p_insert_random, p_delete_random):
    from genomicode import aptamerlib

    # (node, base) -> p_emission
    emission_probs = aptamerlib._calc_emission_probs(
        library, base2emission, p_mismatch)
    # (node1, node2) -> p_transition
    transition_probs = aptamerlib._calc_transition_probs(
        library, p_mismatch, p_insert, p_delete,
        p_insert_random, p_delete_random)

    # Make a list of all nodes.
    all_nodes = {}
    for node1, node2 in transition_probs:
        all_nodes[node1] = 1
        all_nodes[node2] = 1

    # Make a list of all bases.
    all_bases = {}
    for node, base in emission_probs:
        assert node in all_nodes
        all_bases[base] = 1
    all_bases = sorted(all_bases)

    # Sort the nodes.
    all_nodes = sorted(all_nodes)

    # Print the total number of nodes.
    x = 0, "", "META", "TOTAL_NODES", len(all_nodes)
    print "\t".join(map(str, x))

    # Print the number of main nodes.
    x = [x for x in all_nodes if x.node_type == aptamerlib.MAIN]
    x = 0, "", "META", "MAIN_NODES", len(x)
    print "\t".join(map(str, x))

    for i, node1 in enumerate(all_nodes):
        x = "%s-%s-%s-%s" % (
            node1.node_type, node1.i_seqset, node1.i_sequence, node1.i_base)
        x = i+1, node1, "META", "NODE", x
        print "\t".join(map(str, x))
        
        # Count the number of transitions.
        x = [node2 for node2 in all_nodes
             if (node1, node2) in transition_probs]
        x = i+1, node1, "META", "NUM_TRANSITIONS", len(x)
        print "\t".join(map(str, x))

        # Count the unused probability mass.
        p = [transition_probs[(node1, node2)] for node2 in all_nodes
             if (node1, node2) in transition_probs]
        unused = 1.0 - sum(p)
        x = i+1, node1, "META", "UNUSED_P", "%.3g" % unused
        print "\t".join(map(str, x))
        
        for base in all_bases:
            if (node1, base) not in emission_probs:
                continue
            p = emission_probs[(node1, base)]
            p_str = "%.3g" % p
            x = i+1, node1, "EMISSION", base, p_str
            print "\t".join(map(str, x))
        for node2 in all_nodes:
            if (node1, node2) not in transition_probs:
                continue
            p = transition_probs[(node1, node2)]
            p_str = "%.3g" % p
            x = i+1, node1, "TRANSITION", node2, p_str
            print "\t".join(map(str, x))

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
        "--show_markov_model_only", default=False, action="store_true",
        help="Print out the markov model and exit (for debugging).")
    parser.add_argument(
        "--alignment_file", help="File to store the alignments.")
    parser.add_argument(
        "--markov_file", help="File to store the model for each alignment.")
    #parser.add_argument(
    #    "--output_format", choices=[OUT_TABLE, OUT_MARKOV, OUT_ALIGNMENT],
    #    default=OUT_TABLE, 
    #    help="What kind of output to generate.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")

    group = parser.add_argument_group(title="Filter")
    group.add_argument(
        "--max_alignments", type=int,
        help="Maximum number of sequences to align.")
    group.add_argument(
        "--min_seqlen", type=int, default=None,
        help="Discard sequences less than this minimum length.")
    group.add_argument(
        "--titles", default=[], action="append",
        help="Want reads with these titles.  "
        "Comma-separated titles, parameter can be used multiple times.")
    
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
        "--delete", type=float, default=DEF_DELETE,
        help="Probability of a deletion (default %0.2f)." % DEF_DELETE)
    group.add_argument(
        "--insert_random", type=float, default=DEF_INSERT_RANDOM,
        help="Probability of an insertion in the random region "
        "(default %0.2f)." % DEF_INSERT_RANDOM)
    group.add_argument(
        "--delete_random", type=float, default=DEF_DELETE_RANDOM,
        help="Probability of a deletion in the random region "
        "(default %0.2f).  NOT USED." % DEF_DELETE_RANDOM)
    
    args = parser.parse_args()

    assert os.path.exists(args.library_file), \
           "File not found: %s" % args.library_file
    assert os.path.exists(args.sequence_file), \
           "File not found: %s" % args.sequence_file
    assert args.num_procs >= 1 and args.num_procs < 256
    assert args.min_seqlen is None or (
        args.min_seqlen >= 0 and args.min_seqlen < 100), args.min_seqlen

    assert args.mismatch > 0 and args.mismatch < 0.5, \
           "mismatch (%s) should be between 0 and 0.5" % args.mismatch
    assert args.insert > 0 and args.insert < 0.5, \
           "insert (%s) should be between 0 and 0.5" % args.insert
    assert args.delete > 0 and args.delete < 0.5, \
           "delete (%s) should be between 0 and 0.5" % args.delete
    assert args.insert_random > 0 and args.insert_random < 0.5, \
           "insert_random (%s) should be between 0 and 0.5" % \
           args.insert_random
    assert args.delete_random > 0 and args.delete_random < 0.5, \
           "delete_random (%s) should be between 0 and 0.5" % \
           args.delete_random
    titles = _parse_titles(args.titles)

    base2emission = _parse_base2emission(args.lib2seq)

    library = aptamerlib.read_library(args.library_file)
    mm = aptamerlib.make_markov_model(
        library, base2emission, args.mismatch, args.insert, args.delete,
        args.insert_random, args.delete_random)

    if args.show_markov_model_only:
        write_markov_model(
            library, base2emission, args.mismatch, args.insert, args.delete,
            args.insert_random, args.delete_random)
        return

    if args.alignment_file:
        open(args.alignment_file, 'w')
    if args.markov_file:
        open(args.markov_file, 'w')

    manager = multiprocessing.Manager()
    lock = manager.Lock()
    pool = multiprocessing.Pool(args.num_procs)

    header = [
        "Title", "Is Revcomp", "Score", "Barcode", "Random Region",
        "Length", "Total Errors", 
        "Num Matches", "Num Mismatches", "Num Insertions", "Num Deletions",
        "Ideal Sequence", "Observed Sequence"]
    print "\t".join(header)
    sys.stdout.flush()   # needed for multiprocessing
        
    results = []
    for i in range(args.num_procs):
        fn_args = (
            i, args.num_procs, args.sequence_file,
            args.max_alignments, args.min_seqlen, titles,
            mm, library, base2emission,
            args.alignment_file, args.markov_file, lock)
        keywds = {}
        
        if args.num_procs == 1:
            _align_aptamers_h(*fn_args, **keywds)
        else:
            x = pool.apply_async(_align_aptamers_h, fn_args, keywds)
            results.append(x)
    
    pool.close()
    pool.join()
    for x in results:
        x.get()


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
