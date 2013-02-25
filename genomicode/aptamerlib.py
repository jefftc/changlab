"""

Glossary:
sequence   Ordered list of nucleotides, e.g. ACTGCTAA
seqset     Set of alternative sequences.
library    Ordered list of seqsets.

graph            Consists of nodes and edges that indicate topology.
main graph       Nodes and edges for a sequence with no inserts or deletions.
insertion graph  Nodes and edges for insertions.
deletion graph   Nodes and edges for deletions.
 
Each seqset has a name given by the user.  The names are hashed to
guarantee that they are unique.
 
sequence        string
seqset          list of sequences
library         list of seqsets
graph           dictionary of node -> list of next nodes
node            tuple of (node_type, <seqset>, <alternative>, <base>).
<seqset>        Index of seqset.
<alternate>     Index of alternative sequence.
<base>          Index of base of sequence.
START, END      Special nodes for the start and end of the graph.
INSERTEND       Special node that absorbs inserts at end of library.


Functions:
make_markov_model
score_sequence
pretty_sequences

read_library
format_node
parse_node

"""

# _iter_main_graph
# _get_all_nodes
# _get_first_nodes
# _get_next_nodes
# _node_index
# 
# _make_main_graph
# _make_insertion_graph
# _make_deletion_graph
#
# _calc_transition_probs
# _calc_emissions_probs
#
# _add_deletions_to_alignment


MAIN = "MAIN"
INSERT = "INSERT"
START = "START"
END = "END"
INSERTEND = "INSERT-END"

def _iter_main_graph(library):
    for i_seqset, (name, is_random, seqset) in enumerate(library):
        for i_sequence, sequence in enumerate(seqset):
            for i_base, base in enumerate(sequence):
                yield (name, seqset, sequence, base,
                       i_seqset, i_sequence, i_base)


def _get_all_nodes(*graphs):
    # Return a list of all nodes in a graph (multiple graphs).
    nodes = {}
    for graph in graphs:
        for node, next_nodes in graph.iteritems():
            nodes[node] = 1
            for node in next_nodes:
                nodes[node] = 1
    return sorted(nodes)


def _get_first_nodes(library):
    assert len(library)

    nodes = []
    name, is_random, seqset = library[0]
    for i_sequence in range(len(seqset)):
        sequence = seqset[i_sequence]
        assert len(sequence)
        n = MAIN, 0, i_sequence, 0
        nodes.append(n)
    return nodes


def _get_next_nodes(library, i_seqset, i_sequence, i_base):
    assert i_seqset < len(library)
    name, is_random, seqset = library[i_seqset]
    assert i_sequence < len(seqset)
    sequence = seqset[i_sequence]

    next_nodes = []
    # If there are more bases in this sequence, then the next node
    # is the next base in the sequence.
    if i_base < len(sequence)-1:
        n = MAIN, i_seqset, i_sequence, i_base+1
        next_nodes.append(n)
    # If this is the last node of the last seqset, then the next
    # node is the end node.
    elif i_seqset == len(library)-1:
        n = END,
        next_nodes.append(n)
    # If this is the last node of a sequence, then the next node
    # is the first nodes of the next sequence.
    else:
        x, x, next_seqset = library[i_seqset+1]
        for i in range(len(next_seqset)):
            n = MAIN, i_seqset+1, i, 0
            next_nodes.append(n)
    return next_nodes


def _node_index(library, node):
    # Rough measure of the distance of a node from the START node.
    # The index of a node is not unique.  All alternatives in a seqset
    # will have the same indexes.  The distance may be inaccurate if
    # alternatives in a seqset have different lengths.

    # Optimization: Can cache the calculation of iseqset2index for
    # efficiency.
    # index of seqset, START, END, or INSERTEND to the node index.
    iseqset2index = {}

    iseqset2index[START] = 0
    
    node_index = 1
    for i, (name, is_random, seqset) in enumerate(library):
        longest_seq = None
        for seq in seqset:
            if longest_seq is None or len(seq) > longest_seq:
                longest_seq = len(seq)
        assert longest_seq
        iseqset2index[i] = node_index
        node_index += longest_seq
    iseqset2index[END] = node_index
    iseqset2index[INSERTEND] = node_index

    if node[0] in [MAIN, INSERT]:
        i_seqset, i_sequence, i_base = node[1:]
        index = iseqset2index[i_seqset] + i_base
    else:
        assert node[0] in iseqset2index
        index = iseqset2index[node[0]]

    return index
    

def _make_main_graph(library):
    graph = {}  # node -> list of next nodes

    # START points to the first bases of the first seqsets.
    graph[(START,)] = _get_first_nodes(library)

    # Add the other nodes.
    next_nodes = []
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        node = (MAIN, i_seqset, i_sequence, i_base)
        next_nodes = _get_next_nodes(library, i_seqset, i_sequence, i_base)
        assert node not in graph
        graph[node] = next_nodes

    return graph


def _make_insertion_graph(library):
    graph = {}  # node -> list of next nodes

    # START goes to insert to first base of all sequences in the first
    # seqset.
    first_nodes = _get_first_nodes(library)
    insert_nodes = [(INSERT,) + n[1:] for n in first_nodes]
    graph[(START,)] = insert_nodes

    # Inserts point to themselves and the main node.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        main_node = (MAIN, i_seqset, i_sequence, i_base)
        insert_node = (INSERT, i_seqset, i_sequence, i_base)
        assert insert_node not in graph
        graph[insert_node] = [insert_node, main_node]
    # INSERTEND absorbs all remaining bases until the END.
    graph[(INSERTEND,)] = [(INSERTEND,),  (END,)]

    # Each node points to the insert of the next node.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        main_node = (MAIN, i_seqset, i_sequence, i_base)
        next_nodes = _get_next_nodes(library, i_seqset, i_sequence, i_base)
        insert_nodes = []
        for n in next_nodes:
            if n[0] == MAIN:
                nn = (INSERT,) + n[1:]
            elif n[0] == END:
                nn = (INSERTEND,)
            insert_nodes.append(nn)
        assert main_node not in graph
        graph[main_node] = insert_nodes

    return graph
     
    
def _make_deletion_graph(library):
    graph = {}  # state -> list of next states

    # START points to everything.
    nodes = []
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        x = (MAIN, i_seqset, i_sequence, i_base)
        nodes.append(x)
    graph[(START,)] = nodes

    # Make internal edges.
    for x in _iter_main_graph(library):
        x, x, x, x, i_seqset1, i_sequence1, i_base1 = x
        node1 = (MAIN, i_seqset1, i_sequence1, i_base1)
        for x in _iter_main_graph(library):
            x, x, x, x, i_seqset2, i_sequence2, i_base2 = x
            node2 = (MAIN, i_seqset2, i_sequence2, i_base2)
            
            is_next_node = False
            
            if i_seqset1 > i_seqset2:
                pass
            # Each node points to all subsequent nodes in the same
            # sequence.
            elif i_seqset1 == i_seqset2 and \
                i_sequence1 == i_sequence2:
                if i_base1 < i_base2:
                    is_next_node = True
            # Each node points to all nodes of subsequent seqsets.
            elif i_seqset1 < i_seqset2:
                is_next_node = True
                
            if not is_next_node:
                continue
            if node1 not in graph:
                graph[node1] = []
            graph[node1].append(node2)
                
    # Everything points to END.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        node = (MAIN, i_seqset, i_sequence, i_base)
        if node not in graph:
            graph[node] = []
        assert (END,) not in graph[node]
        graph[node].append((END,))

    # Filter out all the main paths (because they're not deletes).
    main_graph = _make_main_graph(library)
    for node, next_nodes in main_graph.iteritems():
        if node not in graph:
            # Missing for last member of the sequences in the last set.
            continue
        x = [x for x in graph[node] if x not in next_nodes]
        graph[node] = x
        if not x:
            del graph[node]

    return graph


def _calc_transition_probs(library, p_main, p_insert, p_delete):
    # Return a dictionary of (node1, node2) -> transition probability.
    import math

    # Transition probabilities.
    assert abs(p_insert + p_delete + p_main - 1.0) < 0.01

    # Make the graph.
    main_graph = _make_main_graph(library)
    insert_graph = _make_insertion_graph(library)
    delete_graph = _make_deletion_graph(library)

    # DEBUG: Print out the topology of the graph.
    #for node in sorted(insert_graph):
    #    next_nodes = sorted(insert_graph[node])
    #    print node, next_nodes
    #import sys; sys.exit(0)

    # Make a list of all nodes in all graphs.
    nodes = _get_all_nodes(main_graph, insert_graph, delete_graph)

    # (start node, next node) -> probability
    probabilities = {}
    for start_node in nodes:
        probs = {}  # next_node -> p
        for next_node in main_graph.get(start_node, []):
            probs[next_node] = p_main
        for next_node in insert_graph.get(start_node, []):
            probs[next_node] = p_insert
        i_start = _node_index(library, start_node)
        for next_node in delete_graph.get(start_node, []):
            i_next = _node_index(library, next_node)
            assert i_next-i_start >= 2
            probs[next_node] = math.pow(p_delete, i_next-i_start-1)
        if not probs:  # no next nodes
            continue

        # Normalize the probabilities to 1.0.
        total = sum(probs.values())
        for next_node in probs:
            probs[next_node] = probs[next_node] / total
        for next_node, p in probs.iteritems():
            probabilities[(start_node, next_node)] = p
    return probabilities

    
def _calc_emission_probs(library, base2emission, p_match, p_mismatch):
    # Make a list of all the emissions.
    emissions = {}
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        e = base2emission.get(base, base)
        emissions[e] = 1
    assert "START" not in emissions
    assert "END" not in emissions
    emissions = ["START", "END"] + sorted(emissions)
    emissions_in_library = emissions[2:]


    probabilities = {}  # (node, emission) -> p
    
    # Set emission probabilities for the main graph.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        node = MAIN, i_seqset, i_sequence, i_base

        probs = {}
        for e in emissions_in_library:
            p = p_match
            if base2emission.get(base, base) != e:
                p = p_mismatch
            probs[e] = p
        # Normalize the probabilities to 1.0.
        total = sum(probs.values())
        for b in probs:
            probs[b] = probs[b] / total
        for b, p in probs.iteritems():
            probabilities[(node, b)] = p
    
    # Set emission probabilities for the insertions.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        node = INSERT, i_seqset, i_sequence, i_base
        for e in emissions_in_library:
            probabilities[(node, e)] = 1.0/len(emissions_in_library)

    # Set emission probabilities for INSERTEND.
    for e in emissions_in_library:
        probabilities[((INSERTEND,), e)] = 1.0/len(emissions_in_library)

    # Set emission probabilities for start and end.
    probabilities[((START,), "START")] = 1.0
    probabilities[((END,), "END")] = 1.0

    return probabilities
             

def make_markov_model(
        library, base2emission, p_insert=None, p_delete=None, p_mismatch=None):
    # Return a MarkovModel.
    import numpy
    from Bio import MarkovModel

    p_mismatch = p_mismatch or 0.10
    p_match = 1.0 - p_mismatch
    assert p_mismatch >= 0 and p_mismatch < 0.50
    assert p_match > 0 and p_match <= 1.0
    
    p_insert = p_insert or p_mismatch
    p_delete = p_delete or (p_mismatch**2)  # transition and mismatch
    p_main = 1.0 - p_insert - p_delete
    assert p_insert >= 0 and p_insert < 0.50
    assert p_delete >= 0 and p_delete < 0.50
    assert p_main >= 0 and p_main <= 1.0

    # Calculate the transition probabilities.
    transition_probs = _calc_transition_probs(
        library, p_main, p_insert, p_delete)
    emission_probs = _calc_emission_probs(
        library, base2emission, p_match, p_mismatch)
    #for (start, end) in sorted(transition_probs):
    #    x = transition_probs[(start, end)], start, end
    #    print "\t".join(map(str, x))
    #import sys; sys.exit(0)

    # Make a list of all the nodes.
    nodes = {}
    for (node1, node2) in transition_probs:
        nodes[node1] = 1
        nodes[node2] = 1
    nodes = sorted(nodes)

    # Make a list of all the emissions.
    emissions = {}
    for (node, emission) in emission_probs:
        emissions[emission] = 1
    emissions = sorted(emissions)

    # N             number of states
    # M             number of emissions
    # p_initial     N-vector of initial starting probabilities.
    # p_emission    N x M matrix.  Each row sums to 1.0.
    # p_transition  N x N matrix.

    # Index the nodes and emissions for the matrices.
    node2i = {}
    for i, s in enumerate(nodes):
        node2i[s] = i
    emission2i = {}
    for i, s in enumerate(emissions):
        emission2i[s] = i
    
    N = len(nodes)
    M = len(emissions)
    p_initial = numpy.zeros(N)
    p_emission = numpy.zeros((N,M))
    p_transition = numpy.zeros((N,N))

    # Start at the START node.
    p_initial[node2i[(START,)]] = 1.0

    # Set the transition probabilities.
    for (start, end) in transition_probs:
        i = node2i[start]
        j = node2i[end]
        p_transition[i][j] = transition_probs[(start, end)]

    # Set the emissions probabilities.
    for (node, emission) in emission_probs:
        i = node2i[node]
        j = emission2i[emission]
        p_emission[i][j] = emission_probs[(node, emission)]
 
    mm = MarkovModel.MarkovModel(
        nodes, emissions, p_initial, p_transition, p_emission)
    return mm


def _add_deletions_to_alignment(library, alignment):
    # Return list of (node, match_type, base, base in sequence).

    # For each seqset, figure out which alternate is used in the
    # alignment.
    iseqset2isequence = {}
    for x in alignment:
        node, match_type, base_in_library, base_in_sequence = x
        if len(node) < 4:
            continue
        i_seqset, i_sequence, i_base = node[1:]
        if i_seqset in iseqset2isequence:
            assert iseqset2isequence[i_seqset] == i_sequence
        iseqset2isequence[i_seqset] = i_sequence
    # If a seqset is completely deleted, then arbitrarily use sequence
    # 0.
    for i in range(len(library)):
        if i not in iseqset2isequence:
            iseqset2isequence[i] = 0

    # Make a list of all the bases in the library.
    all_bases = []
    for i_seqset, (name, is_random, seqset) in enumerate(library):
        i_sequence = iseqset2isequence[i_seqset]
        sequence = seqset[i_sequence]
        for i_base in range(len(sequence)):
            x = i_seqset, i_sequence, i_base
            all_bases.append(x)

    # Make a list of the bases that are matched.
    indexes2align = {}  # Either (i_seqset, i_sequence, i_base) or INSERTEND.
    for align in alignment:
        node = align[0]
        if len(node) >= 4:
            key = node[1:]
        else:
            assert node[0] == INSERTEND
            key = INSERTEND
        if key not in indexes2align:
            indexes2align[key] = []
        indexes2align[key].append(align)

    # Iterate through the library, adding deletions when necessary.
    full_alignment = []
    for i_seqset, (name, is_random, seqset) in enumerate(library):
        i_sequence = iseqset2isequence[i_seqset]
        sequence = seqset[i_sequence]
        for i_base, base in enumerate(sequence):
            x = i_seqset, i_sequence, i_base
            if x in indexes2align:
                full_alignment.extend(indexes2align[x])
            else:
                node = MAIN, i_seqset, i_sequence, i_base
                align = node, "DELETE", base, "-"
                full_alignment.append(align)
    full_alignment.extend(indexes2align.get(INSERTEND, []))
    return full_alignment


def _score_sequence_h(mm, library, base2emission, sequence):
    # Return score, list of (node, match_type, base, base in sequence).
    import math
    from Bio import MarkovModel
    
    # Add "START" and "END" emissions to the ends.
    assert type(sequence) is type("")
    sequence = ["START"] + list(sequence) + ["END"]

    # Score the alignment.
    alignments = MarkovModel.find_states(mm, sequence)
    states, score = alignments[0]
    lscore = math.log(score)

    #print sequence
    #for x in states:
    #    print repr(x)
    #import sys; sys.exit(0)

    # Remove the ["START"] and ["END"] states.
    assert states[0] == (START,)
    assert states[-1] == (END,)
    states = states[1:-1]
    sequence = sequence[1:-1]

    alignment = []
    for i, node in enumerate(states):
        base_in_seq = sequence[i]

        if node[0] == INSERTEND:
            x = node, "INSERT", "-", base_in_seq
            alignment.append(x)
            continue
        
        assert node[0] in [MAIN, INSERT]
        i_seqset, i_sequence, i_base = node[1:]
        name, is_random, seqset = library[i_seqset]
        seq = seqset[i_sequence]
        
        base_in_lib = "-"
        if node[0] == MAIN:
            base_in_lib = seq[i_base]
        match_type = "MATCH"
        if base2emission.get(base_in_lib, base_in_lib) != base_in_seq:
            match_type = "MISMATCH"
        x = node, match_type, base_in_lib, base_in_seq
        alignment.append(x)

    return lscore, alignment


def score_sequence(mm, library, base2emission, sequence):
    # Return score, is_revcomp,
    # list of (node, match_type, base, base in sequence).
    from Bio import Seq

    # Score both the sequence and its reverse complement.  Return the
    # one with the higher score.
    sequence_rc = Seq.Seq(sequence).reverse_complement().tostring()
    x1 = _score_sequence_h(mm, library, base2emission, sequence)
    x2 = _score_sequence_h(mm, library, base2emission, sequence_rc)
    lscore1, alignment1 = x1
    lscore2, alignment2 = x2

    lscore, is_revcomp, alignment = lscore1, False, alignment1
    if lscore2 > lscore1:
        lscore, is_revcomp, alignment = lscore2, True, alignment2

    #for x in alignment:
    #    print repr(x)
    #import sys; sys.exit(0)

    alignment = _add_deletions_to_alignment(library, alignment)
    return lscore, is_revcomp, alignment


def pretty_sequences(alignment):
    ideal_seq = []
    real_seq = []
    prev_i_seqset = None
    for i, x in enumerate(alignment):
        node, match_type, base_in_lib, base_in_seq = x
        
        i_seqset = i_sequence = i_base = None
        if node[0] in [MAIN, INSERT]:
            i_seqset, i_sequence, i_base = node[1:]

        if i_seqset != prev_i_seqset and prev_i_seqset is not None:
            ideal_seq.append("*")
            real_seq.append("*")
        prev_i_seqset = i_seqset
            
        ideal_seq.append(base_in_lib)
        real_seq.append(base_in_seq)
    ideal_seq = "".join(ideal_seq)
    real_seq = "".join(real_seq)
    return ideal_seq, real_seq


def read_library(filename):
    """read_library(filename) -> list of (name, is_random, seqset)"""
    import filelib

    # Format of file is:
    # "Random Region"  "Name"  "Alternate"  ["Alternate", ...]
    # <1/0>  <name>  <alt1>  <alt2>  <altN>
    # <1/0>  <name>  <alt1>  <alt2>  <altN>
    #
    # Columns are separated by tabs.  Each line contains a seqset.  A
    # seqset is represented as a list of sequences.

    handle = filelib.read_cols(filename)
    header = handle.next()

    # Check the header format.
    assert len(header) >= 3, "invalid library format"
    assert header[0].upper() == "Random Region".upper()
    assert header[1].upper() == "Name".upper()
    for i in range(2, len(header)):
        assert header[i].upper() == "Alternate".upper()
    
    library = []  # list of (is_random, name, seqset)
    for cols in handle:
        assert len(cols) >= 3
        is_random = cols[0]
        name = cols[1]
        sequences = [x.strip() for x in cols[2:] if x.strip()]
        assert is_random in ["0", "1"]
        assert len(sequences) >= 1
        x = name, int(is_random), sequences
        library.append(x)
    assert len(library)

    #handle = filelib.read_cols(filename)
    #names = handle.next()
    #matrix = [x for x in handle]
    #
    #for i in range(len(names)):
    #    name = names[i]
    #    sequences = []
    #    for cols in matrix:
    #        if i < len(cols):
    #            sequences.append(cols[i])
    #    sequences = [x for x in sequences if x]
    #    x = name, sequences
    #    library.append(x)

    return library


def format_node(node):
    return ":".join(map(str, node))


def parse_node(node_str):
    x = node_str.split(":")
    assert len(x) >= 1 and len(x) <= 4
    return tuple(x)

