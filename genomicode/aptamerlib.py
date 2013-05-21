"""

Glossary:
sequence         Ordered list of nucleotides, e.g. "ACTGCTAA"
seqset           Set of alternative sequences, e.g. ["AC", "TG", "TT"]
library          Ordered list of (name, is_random, seqset).

graph            Consists of nodes and edges that indicate topology.
main graph       Nodes and edges that model a sequence with no inserts or dels
insertion graph  Nodes and edges that model insertions.
deletion graph   Nodes and edges that model deletions.
 
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
node_type       Constants: START, MAIN, INSERT, END, INSERT-END.


Nodes
-----
(MAIN, <seqset>, <alternative>, <base>)
(INSERT, <seqset>, <alternative>, <base>)
(START,)
(END,)
(INSERTEND,)


Emission probabilities
----------------------
p_mismatch   Non-matching bases.    Default 0.10.

o p_mismatch is the relative probability density that a node in the
  MAIN graph emits a non-matching base.  The rest of the probability
  density is reserved for matches.
o In the INSERT node, each of the bases can be emitted with equal
  probability.
o The START node and the END node emit START and END with 100%
  probability, respectively.

These are technically relative weights, rather than probabilities.
So, with the default weights, a match gets 9 times the density than a
mismatch.  We do it this way because we don't know how many possible
mismatches there are.  So if the alphabet has 3 possible mismatches,
or if it has 10, a match always has 9 times the density.



Transition probabilities
------------------------
p_insert     Inserting one base.    Default (1-p_match).
p_delbase    Deleting one base.     Default (1-p_match)**2 (transition and mm).
p_delchunk   Deleting one chunk.    Default (1-p_match)**2 (transition and mm).

o p_insert is the probability of a transition from any node (including
  an INSERT node) to an INSERT node.
o p_delbase is the probability of a transition from a MAIN node to
  another MAIN node that is one base away.
o p_delchunk is the probability of a transition from a MAIN node to
  another MAIN node that skips an entire seqset.
o The rest of the probability density is reserved for transitions from
  any node to a MAIN node (without deletions).

As with the emissions, these are relative weights, rather than
probabilities.



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
# _node_dist
# 
# _make_main_graph
# _make_insertion_graph
# _make_deletion_graph
#
# _calc_emission_probs
# _calc_transition_probs
#
# _add_deletions_to_alignment


START = "START"
MAIN = "MAIN"
INSERT = "INSERT"
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


def _node_dist(library, start_node, end_node):
    # Return the distance from the start_node to the end_node as a
    # tuple of (total_bases, num_seqset, num_bases).  total_bases is
    # the minimum number of bases from the start_node to the end_node.
    # num_seqset and num_bases is the minimum number of seqsets and
    # bases that need to be traversed.

    # Possible cases:
    # 1.  Both are START or both are END (or INSERTEND).
    # 2.  start_node is START and end_node is END (or INSERTEND).
    # 3.  start_node is START and end_node is internal.
    # 4.  start_node is internal and end_node is END (or INSERTEND).
    # 5.  Nodes are in the same seqset, same sequence.
    # 6.  Nodes are in different seqsets.
    #
    # There is no distance between bases inserted at the end.

    min_bases = []  # minimum number of bases in each seqset.
    for name, is_random, seqset in library:
        seqlens = [len(x) for x in seqset]
        x = min(seqlens)
        min_bases.append(x)

    assert start_node[0] in [MAIN, START, END, INSERT, INSERTEND]
    assert end_node[0] in [MAIN, START, END, INSERT, INSERTEND]

    # Case 1.
    if start_node[0] == START and end_node[0] == START:
        total_bases, num_seqset, num_bases = 0, 0, 0
    elif start_node[0] in ["END", "INSERTEND"] and \
           end_node[0] in ["END", "INSERTEND"]:
        total_bases, num_seqset, num_bases = 0, 0, 0
    # Case 2.
    elif start_node[0] == START and end_node[0] in ["END", "INSERTEND"]:
        total_bases = sum(min_bases)
        num_seqset = len(min_bases)
        num_bases = 0
    # Case 3.
    elif start_node[0] == START:
        i_seqset2, i_sequence2, i_base2 = end_node[1:]
        total_bases = sum(min_bases[:i_seqset2]) + i_base2
        num_seqset = i_seqset2
        num_bases = i_base2
    # Case 4.
    elif end_node[0] in ["END", "INSERTEND"]:
        i_seqset1, i_sequence1, i_base1 = start_node[1:]
        sequence1 = library[i_seqset1][-1][i_sequence1]
        num_seqset = len(library) - i_seqset1
        num_bases = len(sequence1) - i_base1
        total_bases = sum(min_bases[(i_seqset1+1):]) + num_bases
    # Case 5:
    elif start_node[1] == end_node[1]:
        i_seqset1, i_sequence1, i_base1 = start_node[1:]
        i_seqset2, i_sequence2, i_base2 = end_node[1:]
        assert i_sequence1 == i_sequence2
        assert i_base2 >= i_base1
        num_bases = i_base2 - i_base1
        num_seqset = 0
        total_bases = num_bases
    # Case 6.
    else:
        i_seqset1, i_sequence1, i_base1 = start_node[1:]
        i_seqset2, i_sequence2, i_base2 = end_node[1:]
        assert i_seqset2 > i_seqset1
        sequence1 = library[i_seqset1][-1][i_sequence1]
        num_seqset = i_seqset2 - i_seqset1
        num_bases = (len(sequence1) - i_base1) + i_base2
        total_bases = sum(min_bases[(i_seqset1+1):i_seqset2]) + num_bases

    return total_bases, num_seqset, num_bases


## def _node_index(library, node):
##     # Rough measure of the distance of a node from the START node.
##     # The index of a node is not unique.  All alternatives in a seqset
##     # will have the same indexes.  The distance may be inaccurate if
##     # alternatives in a seqset have different lengths.

##     # Optimization: Can cache the calculation of iseqset2index for
##     # efficiency.
##     # index of seqset, START, END, or INSERTEND to the node index.
##     iseqset2index = {}

##     iseqset2index[START] = 0
    
##     node_index = 1
##     for i, (name, is_random, seqset) in enumerate(library):
##         longest_seq = None
##         for seq in seqset:
##             if longest_seq is None or len(seq) > longest_seq:
##                 longest_seq = len(seq)
##         assert longest_seq
##         iseqset2index[i] = node_index
##         node_index += longest_seq
##     iseqset2index[END] = node_index
##     iseqset2index[INSERTEND] = node_index

##     if node[0] in [MAIN, INSERT]:
##         i_seqset, i_sequence, i_base = node[1:]
##         index = iseqset2index[i_seqset] + i_base
##     else:
##         assert node[0] in iseqset2index
##         index = iseqset2index[node[0]]

##     return index
    

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
    graph[(INSERTEND,)] = [(INSERTEND,), (END,)]

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

    # START points to everything.  Any number of the initial bases can
    # be deleted.
    nodes = []
    for x in _iter_main_graph(library):
        x, x, x, x, i_seqset, i_sequence, i_base = x
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
        x, x, x, x, i_seqset, i_sequence, i_base = x
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


def _calc_emission_probs(library, base2emission, p_mismatch):
    # Make a list of all the bases that can be emitted.
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
        emission = base2emission.get(base, base)  # what this base emits

        probs = {}   # base -> probability
        for e in emissions_in_library:
            p = 1.0-p_mismatch
            if e != emission:
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
             

def _calc_transition_probs(library, p_insert, p_delbase, p_delchunk):
    # Return a dictionary of (node1, node2) -> transition probability.
    import math

    # Transition probabilities.
    p_main = 1.0 - (p_insert + p_delbase + p_delchunk)
    assert p_insert >= 0 and p_insert < 0.5
    assert p_delbase >= 0 and p_delbase < 0.5
    assert p_delchunk >= 0 and p_delchunk < 0.5
    assert p_main > 0 and p_main <= 1.0

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
            assert next_node[0] in ["MAIN", "END", "INSERT", "INSERT-END"], \
                   repr(next_node)
            p = p_main
            if next_node[0] in ["INSERT", "INSERT-END"]:
                p = p_insert
            probs[next_node] = p
        for next_node in delete_graph.get(start_node, []):
            x = _node_dist(library, start_node, next_node)
            total_bases, num_seqsets, num_bases = x
            # If deleting a seqset is lower probability than deleting
            # the individual bases, then just model as deleting the
            # individual bases.
            p1 = pow(p_delbase, total_bases)
            p2 = (pow(p_delbase, num_bases) *
                  pow(p_delchunk, num_seqsets))
            probs[next_node] = max(p1, p2)
        if not probs:  # no next nodes
            continue
        #print "MAIN", main_graph.get(start_node, [])
        #print "INSERT", insert_graph.get(start_node, [])
        #print "DELETE", delete_graph.get(start_node, [])
        #print start_node, probs

        # If there are any non-insertion or deletion nodes, distribute
        # remaining probability mass to them.
        ## num_nodes = len([x for x in probs.itervalues() if x is None])
        ## if num_nodes:
        ##     p = sum([x for x in probs.itervalues() if x is not None])
        ##     p_main = (1.0 - p) / num_nodes
        ##     for next_node in probs:
        ##         if probs[next_node] is None:
        ##             probs[next_node] = p_main
        ##     #print num_nodes, p_main
        ##     #if p_main < 0.20:
        ##     #    for n, x in probs.iteritems():
        ##     #        print n, x

        ## # Otherwise, normalize the probabilities to 1.0.  This can
        ## # happen, for example, if we decide to treat transitions into
        ## # and out of INSERT nodes as insertions.
        ## else:
        # Normalize the probabilities to 1.0.
        total = sum(probs.values())
        for next_node in probs:
            probs[next_node] = probs[next_node] / total
        #print start_node, probs
        #import sys; sys.exit(0)
        #if len(main_graph.get(start_node, [])) >= 4:
        #    print "HERE"
        #    for n, x in probs.iteritems():
        #        print n, x
        #    print "DONE"

        # Make sure probabilities sum to 1.
        #assert abs(sum(probs.values()) - 1.0) < 1E-10
        
        for next_node, p in probs.iteritems():
            probabilities[(start_node, next_node)] = p
    return probabilities

    
def make_markov_model(
    library, base2emission, p_mismatch=None, p_insert=None,
    p_delbase=None, p_delchunk=None):
    # Return a MarkovModel.
    import numpy
    from genomicode import MarkovModel

    if p_mismatch is None:
        p_mismatch = 0.1
    if p_insert is None:
        p_insert = p_mismatch
    if p_delbase is None:
        p_delbase = p_mismatch**2   # transition and mismatch
    if p_delchunk is None:
        p_delchunk = p_mismatch**2  # transition and mismatch
    assert p_mismatch >= 0 and p_mismatch < 0.50
    assert p_insert >= 0 and p_insert < 0.50
    assert p_delbase >= 0 and p_delbase < 0.50
    assert p_delchunk >= 0 and p_delchunk < 0.50

    # Calculate the transition probabilities.
    transition_probs = _calc_transition_probs(
        library, p_insert, p_delbase, p_delchunk)
    emission_probs = _calc_emission_probs(library, base2emission, p_mismatch)
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
    #p_initial = [0.0] * N
    #p_emission = [[0.0]*M for i in range(N)]
    #p_transition = [[0.0]*N for i in range(N)]
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
    from genomicode import MarkovModel
    
    # Add "START" and "END" emissions to the ends.
    assert type(sequence) is type("")
    sequence = ["START"] + list(sequence) + ["END"]

    # Score the alignment.
    alignments = MarkovModel.find_states(mm, sequence)
    states, score = alignments[0]
    score = max(score, 1E-100)
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


def guess_sequence_orientation(sequence, library):
    # Return 1, 0, or -1.  1 if the sequence is in the right
    # orientation, -1 if it needs to be reverse complemented, and 0 if
    # I can't tell.
    from Bio import Seq
    from Bio import pairwise2

    sequence_rc = Seq.Seq(sequence).reverse_complement().tostring()

    # First, see if it matches the first seqset of the library exactly.
    for seq in library[0][-1]:
        if sequence[:len(seq)] == seq:
            return 1
        if sequence_rc[:len(seq)] == seq:
            return -1

    # If it doesn't match exactly, see if it matches to within 2
    # mismatches.
    for seq in library[0][-1]:
        if len(seq) < 6:  # if sequence too short, hard to align.
            continue
        score = pairwise2.align.globalxx(
            sequence[:len(seq)], seq, score_only=True)
        if score >= len(seq)-2:
            return 1
        score = pairwise2.align.globalxx(
            sequence_rc[:len(seq)], seq, score_only=True)
        if score >= len(seq)-2:
            return -1

    # I can't tell.
    return 0


def score_sequence(mm, library, base2emission, sequence):
    # Return score, is_revcomp,
    # list of (node, match_type, base, base in sequence).
    from Bio import Seq
    from Bio import pairwise2

    assert library

    orientation = guess_sequence_orientation(sequence, library)
    
    if orientation == 1:
        x = _score_sequence_h(mm, library, base2emission, sequence)
        lscore, alignment = x
        is_revcomp = False
    elif orientation == -1:
        sequence_rc = Seq.Seq(sequence).reverse_complement().tostring()
        x = _score_sequence_h(mm, library, base2emission, sequence_rc)
        lscore, alignment = x
        is_revcomp = True
    else:
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


