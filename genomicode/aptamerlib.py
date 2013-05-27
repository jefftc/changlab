"""

Glossary:
sequence         Ordered list of nucleotides, e.g. "ACTGCTAA"
seqset           Set of alternative sequences, e.g. ["AC", "TG", "TT"]
library          Ordered list of seqsets.

Each seqset has a name and other data (whether it is a barcode,
whether it is part of the random region) given by the user.  The names
are hashed to guarantee that they are unique.
 
graph            Consists of nodes and edges that indicate topology.
main graph       Nodes and edges that model a sequence with no inserts or dels
insertion graph  Nodes and edges that model insertions.
deletion graph   Nodes and edges that model deletions.

sequence         string
seqset           SequenceSet object
library          list of seqsets
graph            dictionary of node -> list of next nodes
node             Node object with node_type, <seqset>, <alternate>, <base>
node_type        Constants: START, MAIN, INSERT, END, INSERTEND.
<seqset>         Index of seqset.
<alternate>      Index of alternative sequence.
<base>           Index of base of sequence.
START, END       Special nodes for the start and end of the graph.
INSERTEND        Special node that absorbs inserts at end of library.


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
p_delbase    Deleting one base.     Default (1-p_match)*p_match*p_match.
p_delchunk   Deleting one chunk.    Default p_delbase * p_match*p_match.

o p_insert is the probability of a transition from any node (including
  an INSERT node) to an INSERT node.
o p_delbase is the probability of a transition from a MAIN node to
  another MAIN node that is one base away.
  Default should be p_trans * p_mismatch * p_trans.
o p_delchunk is the probability of a transition from a MAIN node to
  another MAIN node that skips an entire seqset.
  Default should be p_del * p_match * p_trans (chunk with one base
  deletion, and one base match).
o The rest of the probability density is reserved for transitions from
  any node to a MAIN node (without deletions).

As with the emissions, these are relative weights, rather than
probabilities.



Functions:
make_markov_model
score_sequence
pretty_sequences

parse_node

read_library
parse_fastq

"""

class Node:
    def __init__(self, node_type, i_seqset=None, i_sequence=None, i_base=None):
        assert node_type in [MAIN, START, END, INSERT, INSERTEND]
        self.node_type = node_type
        self.i_seqset = i_seqset
        self.i_sequence = i_sequence
        self.i_base = i_base
    def __str__(self):
        x = [self.node_type]
        if self.i_seqset is not None:
            x += [self.i_seqset, self.i_sequence, self.i_base]
        return "Node(%s)" % (", ".join(map(str, x)))
    def __hash__(self):
        return hash(
            (self.node_type, self.i_seqset, self.i_sequence, self.i_base))
    def __eq__(self, other):
        return (self.node_type == other.node_type and 
                self.i_seqset == other.i_seqset and 
                self.i_sequence == other.i_sequence and 
                self.i_base == other.i_base)
    def __cmp__(self, other):
        x1 = self.node_type, self.i_seqset, self.i_sequence, self.i_base
        x2 = other.node_type, other.i_seqset, other.i_sequence, other.i_base
        return cmp(x1, x2)



class SequenceSet:
    def __init__(self, name, is_barcode, is_random, alternates):
        self.name = name
        self.is_barcode = is_barcode
        self.is_random = is_random
        self.alternates = alternates[:]


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
    for i_seqset, seqset in enumerate(library):
        name = library[i_seqset].name
        for i_sequence, sequence in enumerate(seqset.alternates):
            for i_base, base in enumerate(sequence):
                yield (name, seqset, sequence, base,
                       i_seqset, i_sequence, i_base)


def _node2sortkey(node):
    # (order1, i_seqset, i_sequence, i_base, order2)
    nt2order1 = {
        START : 0,
        INSERT : 1,
        MAIN : 1,
        INSERTEND : 2,
        END : 3,
        }
    nt2order2 = {
        START : 0,
        INSERT : 0,
        MAIN : 1,
        INSERTEND : 0,
        END : 0,
        }
    assert node.node_type in nt2order1
    assert node.node_type in nt2order2
    order1 = nt2order1[node.node_type]
    order2 = nt2order2[node.node_type]
    x = order1, node.i_seqset, node.i_sequence, node.i_base, order2
    return x
    

def _get_all_nodes(*graphs):
    # Return a list of all nodes in a graph (multiple graphs).
    nodes = {}
    for graph in graphs:
        for node, next_nodes in graph.iteritems():
            nodes[node] = 1
            for node in next_nodes:
                nodes[node] = 1
    return sorted(nodes, key=_node2sortkey)


def _get_first_nodes(library):
    assert len(library)

    nodes = []
    seqset = library[0]
    for i_sequence, sequence in enumerate(seqset.alternates):
        assert len(sequence)
        n = Node(MAIN, 0, i_sequence, 0)
        nodes.append(n)
    return nodes


def _get_next_nodes(library, i_seqset, i_sequence, i_base):
    assert i_seqset < len(library)
    seqset = library[i_seqset]
    assert i_sequence < len(seqset.alternates)
    sequence = seqset.alternates[i_sequence]

    next_nodes = []
    # If there are more bases in this sequence, then the next node
    # is the next base in the sequence.
    if i_base < len(sequence)-1:
        n = Node(MAIN, i_seqset, i_sequence, i_base+1)
        next_nodes.append(n)
    # If this is the last node of the last seqset, then the next
    # node is the end node.
    elif i_seqset == len(library)-1:
        n = Node(END)
        next_nodes.append(n)
    # If this is the last node of a sequence, then the next node
    # is the first nodes of the next sequence.
    else:
        next_seqset = library[i_seqset+1]
        for i in range(len(next_seqset.alternates)):
            n = Node(MAIN, i_seqset+1, i, 0)
            next_nodes.append(n)
    return next_nodes


def _node_dist(library, start_node, end_node):
    # Return the distance from the start_node to the end_node as a
    # tuple of (total_bases, num_seqset, num_bases).  total_bases is
    # the minimum number of bases from the start_node to the end_node.
    # num_seqset and num_bases is the minimum number of (random
    # region) seqsets and bases that need to be traversed.

    # Possible cases:
    # 1.  Both are START or both are END (or INSERTEND).
    # 2.  Nodes are in the same seqset, same sequence.
    # 3.  start_node is START and end_node is END (or INSERTEND).
    # 4.  start_node is START and end_node is internal.
    # 5.  start_node is internal and end_node is END (or INSERTEND).
    # 6.  start_node and end_node are internal.
    #
    # There is no distance between bases inserted at the end.

    min_bases = []  # minimum number of bases in each seqset.
    for seqset in library:
        seqlens = [len(x) for x in seqset.alternates]
        x = min(seqlens)
        min_bases.append(x)

    assert start_node.node_type in [MAIN, START, END, INSERT, INSERTEND]
    assert end_node.node_type in [MAIN, START, END, INSERT, INSERTEND]

    # Case 1.
    if start_node.node_type == START and end_node.node_type == START:
        total_bases, num_seqset, num_bases = 0, 0, 0
    elif start_node.node_type in ["END", "INSERTEND"] and \
           end_node.node_type in ["END", "INSERTEND"]:
        total_bases, num_seqset, num_bases = 0, 0, 0
    # Case 2:
    elif start_node.i_seqset == end_node.i_seqset:
        assert start_node.i_sequence == end_node.i_sequence
        assert start_node.i_base < end_node.i_base
        num_bases = end_node.i_base - start_node.i_base
        num_seqset = 0
        total_bases = num_bases
    # Case 3-6.
    else:
        total_bases = num_bases = num_seqset = 0
        start_i = 0
        if start_node.node_type != START:
            start_i = start_node.i_seqset+1
            seqset = library[start_node.i_seqset]
            seq = seqset.alternates[start_node.i_sequence]
            num_bases += len(seq)-start_node.i_base
            total_bases += len(seq)-start_node.i_base
        end_i = len(library)
        if end_node.node_type != END:
            end_i = end_node.i_seqset
            num_bases += end_node.i_base + 1   # base 0 is 1 base
            total_bases += end_node.i_base + 1
        for i in range(start_i, end_i):
            seqset = library[i]
            if seqset.is_random:
                num_seqset += 1
            else:
                num_bases += min_bases[i]
            total_bases += min_bases[i]
            
    return total_bases, num_seqset, num_bases


def _make_main_graph(library):
    graph = {}  # node (as string) -> list of next nodes

    # START points to the first bases of the first seqsets.
    graph[Node(START)] = _get_first_nodes(library)

    # Add the other nodes.
    next_nodes = []
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        node = Node(MAIN, i_seqset, i_sequence, i_base)
        next_nodes = _get_next_nodes(library, i_seqset, i_sequence, i_base)
        assert node not in graph
        graph[node] = next_nodes

    return graph


def _make_insertion_graph(library):
    graph = {}  # node (as string) -> list of next nodes

    # START goes to insert to first base of all sequences in the first
    # seqset.
    first_nodes = _get_first_nodes(library)
    insert_nodes = [
        Node(INSERT, n.i_seqset, n.i_sequence, n.i_base) for n in first_nodes]
    graph[Node(START)] = insert_nodes

    # Inserts point to themselves and the main node.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        main_node = Node(MAIN, i_seqset, i_sequence, i_base)
        insert_node = Node(INSERT, i_seqset, i_sequence, i_base)
        assert insert_node not in graph
        graph[insert_node] = [insert_node, main_node]
    # INSERTEND absorbs all remaining bases until the END.
    graph[Node(INSERTEND)] = [Node(INSERTEND), Node(END)]

    # Each node points to the insert of the next node.
    for x in _iter_main_graph(library):
        name, seqset, sequence, base, i_seqset, i_sequence, i_base = x
        main_node = Node(MAIN, i_seqset, i_sequence, i_base)
        next_nodes = _get_next_nodes(library, i_seqset, i_sequence, i_base)
        insert_nodes = []
        for n in next_nodes:
            if n.node_type == MAIN:
                nn = Node(INSERT, n.i_seqset, n.i_sequence, n.i_base)
            elif n.node_type == END:
                nn = Node(INSERTEND, n.i_seqset, n.i_sequence, n.i_base)
            insert_nodes.append(nn)
        assert main_node not in graph
        graph[main_node] = insert_nodes

    return graph
     
    
def _make_deletion_graph(library):
    graph = {}  # node -> list of next nodes

    # START points to everything.  Any number of the initial bases can
    # be deleted.
    nodes = []
    for x in _iter_main_graph(library):
        x, x, x, x, i_seqset, i_sequence, i_base = x
        x = Node(MAIN, i_seqset, i_sequence, i_base)
        nodes.append(x)
    graph[Node(START)] = nodes

    # Make internal edges.
    for x in _iter_main_graph(library):
        x, x, x, x, i_seqset1, i_sequence1, i_base1 = x
        node1 = Node(MAIN, i_seqset1, i_sequence1, i_base1)
        for x in _iter_main_graph(library):
            x, x, x, x, i_seqset2, i_sequence2, i_base2 = x
            node2 = Node(MAIN, i_seqset2, i_sequence2, i_base2)
            
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
        node = Node(MAIN, i_seqset, i_sequence, i_base)
        if node not in graph:
            graph[node] = []
        assert Node(END) not in graph[node]
        graph[node].append(Node(END))

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
    # Can cache _iter_main_graph to optimize.
    
    # Make a list of all the bases that can be emitted.
    emissions = {}
    for x in _iter_main_graph(library):
        x, x, x, base, x, x, x = x
        e = base2emission.get(base, base)
        e = e.upper()   # do case insensitive search
        emissions[e] = 1
    assert "START" not in emissions
    assert "END" not in emissions
    emissions_in_library = sorted(emissions)
    emissions = ["START", "END"] + emissions_in_library

    probabilities = {}  # (node, emission) -> p
    
    # Set emission probabilities for the main graph.
    for x in _iter_main_graph(library):
        x, x, x, base, i_seqset, i_sequence, i_base = x
        node = Node(MAIN, i_seqset, i_sequence, i_base)
        emission = base2emission.get(base, base)  # what this base emits
        emission = emission.upper()

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
        x, x, x, x, i_seqset, i_sequence, i_base = x
        node = Node(INSERT, i_seqset, i_sequence, i_base)
        for e in emissions_in_library:
            probabilities[(node, e)] = 1.0/len(emissions_in_library)

    # Set emission probabilities for INSERTEND.
    for e in emissions_in_library:
        probabilities[(Node(INSERTEND), e)] = 1.0/len(emissions_in_library)

    # Set emission probabilities for start and end.
    probabilities[(Node(START), "START")] = 1.0
    probabilities[(Node(END), "END")] = 1.0

    return probabilities
             

def _calc_transition_probs(library, p_insert, p_delbase, p_delchunk):
    # Return a dictionary of (node1, node2) -> transition probability.

    # Transition probabilities.
    assert p_insert >= 0 and p_insert < 0.5
    assert p_delbase >= 0 and p_delbase < 0.5
    assert p_delchunk >= 0 and p_delchunk < 0.5

    # Make the graph.
    main_graph = _make_main_graph(library)
    insert_graph = _make_insertion_graph(library)
    delete_graph = _make_deletion_graph(library)

    # Make a list of all nodes in all graphs.
    nodes = _get_all_nodes(main_graph, insert_graph, delete_graph)

    # Calculate all the probabilities.  The probabilities will not be
    # normalized and may not sum to 1.
    probabilities = {}   # node -> next_node -> probability.
    for start_node in nodes:
        p_main = 1.0 - (p_insert + p_delbase + p_delchunk)
        if start_node.node_type == INSERT:
            # INSERT nodes do not go to DELETE.
            p_main = 1.0 - p_insert

        probs = {}  # next_node -> p
        for next_node in main_graph.get(start_node, []):
            probs[next_node] = p_main
        for next_node in insert_graph.get(start_node, []):
            assert next_node.node_type in [MAIN, END, INSERT, INSERTEND]
            p = p_main
            if next_node.node_type in [INSERT, INSERTEND]:
                p = p_insert
            probs[next_node] = p
        for next_node in delete_graph.get(start_node, []):
            assert next_node not in probs
            x = _node_dist(library, start_node, next_node)
            total_bases, num_seqsets, num_bases = x
            #print start_node, next_node, total_bases, num_seqsets, num_bases
            assert total_bases >= 1
            assert num_bases >= 1
            # If deleting a seqset is lower probability than deleting
            # the individual bases, then just model as deleting the
            # individual bases.
            p1 = pow(p_delbase, total_bases-1)
            p2 = (pow(p_delbase, num_bases-1) *
                  pow(p_delchunk, num_seqsets))
            probs[next_node] = max(p1, p2)
        if not probs:  # no next nodes
            continue
        probabilities[start_node] = probs


    # Normalize the probabilities.

    # If we normalize the probabilities for each node, then later
    # deletions will have higher probability than earlier deletions
    # (because they have fewer transitions).  The algorithm will favor
    # deleting later bases.
    #
    # To solve this, keep track of the normalization for the first
    # node in the random region.  This should have the most
    # transitions.  Normalize all the nodes in the random region the
    # same way.

    # The START node has the most transitions.  Transitions:
    # 1.  To all MAIN nodes.
    # 2.  To INSERT of MAIN node in seqset 0.
    # 3.  To INSERTEND (not currently used).
    # 4.  To END (not currently used).

    normalized = {}  # start_node -> next_nodes -> p
    for start_node in probabilities:
        probs = probabilities[start_node]
        total = sum(probs.values())
        for next_node in probs:
            probs[next_node] = probs[next_node] / total

        # Make sure probabilities sum to 1.
        assert abs(sum(probs.values()) - 1.0) < 1E-10
        # Make sure probabilities sum to <= 1.
        #total = sum(probs.values())
        #assert total <= 1+1E-8, "%s Prob: %s" % (start_node, total)
            
        normalized[start_node] = probs

    clean = {}
    for start_node in probabilities:
        probs = probabilities[start_node]
        for next_node, p in probs.iteritems():
            clean[(start_node, next_node)] = p
    return clean
    
        
    ##     # Normalize the probabilities.  If this is not in the random
    ##     # region, normalize the probabilities so that they sum to 1.
    ##     # If it is the random region, normalize the probabilities for
    ##     # all nodes the same way.
    ##     is_random = False
    ##     if start_node.node_type == MAIN and \
    ##            library[start_node.i_seqset].is_random:
    ##         is_random = True
    ##     if not is_random or start_node.node_type == INSERT:
    ##         total = sum(probs.values())
    ##         nf = 1.0 / total
    ##     elif norm_factor is None:
    ##         total = sum(probs.values())
    ##         norm_factor = 1.0 / total
    ##         nf = norm_factor
    ##     else:
    ##         nf = norm_factor

    ##     if is_random:
    ##         nns = sorted(probs, key=_node2sortkey)
    ##         for n in nns:
    ##             x = "HERE", start_node, n, probs[n]
    ##             print "\t".join(map(str, x))

    ##     for next_node in probs:
    ##         probs[next_node] = probs[next_node] * nf

    ##     # Make sure probabilities sum to 1.
    ##     #assert abs(sum(probs.values()) - 1.0) < 1E-10
    ##     # Make sure probabilities sum to <= 1.
    ##     total = sum(probs.values())
    ##     assert total <= 1+1E-8, "%s Prob: %s" % (start_node, total)
        
    ##     for next_node, p in probs.iteritems():
    ##         probabilities[(start_node, next_node)] = p


    ## # Make sure START is the first node.
    ## #assert nodes
    ## #assert nodes[0].node_type == START
    
    ## # (start node, next node) -> probability
    ## norm_factor = None
    ## probabilities = {}
    ## for start_node in nodes:
    ##     p_main = 1.0 - (p_insert + p_delbase + p_delchunk)
    ##     if start_node.node_type == INSERT:
    ##         # INSERT nodes do not go to DELETE.
    ##         p_main = 1.0 - p_insert

    ##     probs = {}  # next_node -> p
    ##     for next_node in main_graph.get(start_node, []):
    ##         probs[next_node] = p_main
    ##     for next_node in insert_graph.get(start_node, []):
    ##         assert next_node.node_type in [MAIN, END, INSERT, INSERTEND]
    ##         p = p_main
    ##         if next_node.node_type in [INSERT, INSERTEND]:
    ##             p = p_insert
    ##         probs[next_node] = p

    ##     for next_node in delete_graph.get(start_node, []):
    ##         assert next_node not in probs
    ##         x = _node_dist(library, start_node, next_node)
    ##         total_bases, num_seqsets, num_bases = x
    ##         #print start_node, next_node, total_bases, num_seqsets, num_bases
    ##         assert total_bases >= 1
    ##         assert num_bases >= 1
    ##         # If deleting a seqset is lower probability than deleting
    ##         # the individual bases, then just model as deleting the
    ##         # individual bases.
    ##         p1 = pow(p_delbase, total_bases-1)
    ##         p2 = (pow(p_delbase, num_bases-1) *
    ##               pow(p_delchunk, num_seqsets))
    ##         probs[next_node] = max(p1, p2)
    ##     if not probs:  # no next nodes
    ##         continue

    ##     # Normalize the probabilities.  If this is not in the random
    ##     # region, normalize the probabilities so that they sum to 1.
    ##     # If it is the random region, normalize the probabilities for
    ##     # all nodes the same way.
    ##     is_random = False
    ##     if start_node.node_type == MAIN and \
    ##            library[start_node.i_seqset].is_random:
    ##         is_random = True
    ##     if not is_random or start_node.node_type == INSERT:
    ##         total = sum(probs.values())
    ##         nf = 1.0 / total
    ##     elif norm_factor is None:
    ##         total = sum(probs.values())
    ##         norm_factor = 1.0 / total
    ##         nf = norm_factor
    ##     else:
    ##         nf = norm_factor

    ##     if is_random:
    ##         nns = sorted(probs, key=_node2sortkey)
    ##         for n in nns:
    ##             x = "HERE", start_node, n, probs[n]
    ##             print "\t".join(map(str, x))

    ##     for next_node in probs:
    ##         probs[next_node] = probs[next_node] * nf

    ##     # Make sure probabilities sum to 1.
    ##     #assert abs(sum(probs.values()) - 1.0) < 1E-10
    ##     # Make sure probabilities sum to <= 1.
    ##     total = sum(probs.values())
    ##     assert total <= 1+1E-8, "%s Prob: %s" % (start_node, total)
        
    ##     for next_node, p in probs.iteritems():
    ##         probabilities[(start_node, next_node)] = p
    ## return probabilities

    
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
    p_initial[node2i[Node(START)]] = 1.0

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
        #print repr(x)
        node, log_score, match_type, base_in_library, base_in_sequence = x
        if node.node_type in [START, END, INSERTEND]:
            continue
        if node.i_seqset in iseqset2isequence:
            assert iseqset2isequence[node.i_seqset] == node.i_sequence
        iseqset2isequence[node.i_seqset] = node.i_sequence
    # If a seqset is completely deleted, then arbitrarily use sequence
    # 0.
    for i in range(len(library)):
        if i not in iseqset2isequence:
            iseqset2isequence[i] = 0

    # Make a list of all the bases in the library.
    all_bases = []
    for i_seqset, seqset in enumerate(library):
        i_sequence = iseqset2isequence[i_seqset]
        sequence = seqset.alternates[i_sequence]
        for i_base in range(len(sequence)):
            x = i_seqset, i_sequence, i_base
            all_bases.append(x)

    # Make a list of the bases that are matched.
    indexes2align = {}  # Either (i_seqset, i_sequence, i_base) or INSERTEND.
    for align in alignment:
        node = align[0]
        if node.node_type in [MAIN, INSERT]:
            key = node.i_seqset, node.i_sequence, node.i_base
        else:
            assert node.node_type == INSERTEND
            key = INSERTEND
        if key not in indexes2align:
            indexes2align[key] = []
        indexes2align[key].append(align)

    # Iterate through the library, adding deletions when necessary.
    full_alignment = []
    for i_seqset, seqset in enumerate(library):
        i_sequence = iseqset2isequence[i_seqset]
        sequence = seqset.alternates[i_sequence]
        for i_base, base in enumerate(sequence):
            x = i_seqset, i_sequence, i_base
            if x in indexes2align:
                full_alignment.extend(indexes2align[x])
            else:
                node = Node(MAIN, i_seqset, i_sequence, i_base)
                align = node, 0, "DELETE", base, "-"
                full_alignment.append(align)
    full_alignment.extend(indexes2align.get(INSERTEND, []))
    return full_alignment


def _score_sequence_h(mm, library, base2emission, sequence):
    # Return log(score), list of
    #   (node, log(score), match_type, base in lib, base in sequence).
    import math
    from genomicode import MarkovModel
    
    # Add "START" and "END" emissions to the ends.
    assert type(sequence) is type("")
    sequence = ["START"] + list(sequence) + ["END"]
    sequence_u = [x.upper() for x in sequence]   # score case insensitive

    # Score the alignment.
    alignments = MarkovModel.find_states(mm, sequence_u)
    states, scores = alignments[0]
    log_scores = [math.log(max(x, 1E-100)) for x in scores]
    lscore = log_scores[-1]

    #print sequence
    #for i in range(len(states)):
    #    print states[i], scores[i]
    #import sys; sys.exit(0)

    # Remove the ["START"] and ["END"] states.
    assert states[0] == Node(START)
    assert states[-1] == Node(END)
    states = states[1:-1]
    scores = scores[1:-1]
    log_scores = log_scores[1:-1]
    sequence = sequence[1:-1]

    alignment = []
    for i, node in enumerate(states):
        base_in_seq = sequence[i]

        if node.node_type == INSERTEND:
            x = node, lscore, "INSERT", "-", base_in_seq
            alignment.append(x)
            continue
        
        assert node.node_type in [MAIN, INSERT]
        seq = library[node.i_seqset].alternates[node.i_sequence]
        
        base_in_lib = "-"
        if node.node_type == MAIN:
            base_in_lib = seq[node.i_base]
        b1, b2 = base2emission.get(base_in_lib, base_in_lib), base_in_seq
        match_type = "MATCH"
        if node.node_type == "INSERT":
            match_type = "INSERT"
        elif b1.upper() != b2.upper():
            match_type = "MISMATCH"
        x = node, log_scores[i], match_type, base_in_lib, base_in_seq
        alignment.append(x)

    return lscore, alignment


def guess_sequence_orientation(sequence, library):
    # Return 1, 0, or -1.  1 if the sequence is in the right
    # orientation, -1 if it needs to be reverse complemented, and 0 if
    # I can't tell.
    from Bio import Seq
    from Bio import pairwise2

    sequence = sequence.upper()
    sequence_rc = Seq.Seq(sequence).reverse_complement().tostring()

    # First, see if it matches the first seqset of the library exactly.
    for seq in library[0].alternates:
        if sequence[:len(seq)] == seq.upper():
            return 1
        if sequence_rc[:len(seq)] == seq.upper():
            return -1

    # If it doesn't match exactly, see if it matches to within 2
    # mismatches.
    for seq in library[0].alternates:
        if len(seq) < 6:  # if sequence too short, hard to align.
            continue
        score = pairwise2.align.globalxx(
            sequence[:len(seq)], seq.upper(), score_only=True)
        if score >= len(seq)-2:
            return 1
        score = pairwise2.align.globalxx(
            sequence_rc[:len(seq)], seq.upper(), score_only=True)
        if score >= len(seq)-2:
            return -1

    # I can't tell.
    return 0


def score_sequence(mm, library, base2emission, sequence):
    # Return score, is_revcomp,
    # list of (node, match_type, base, base in sequence).
    from Bio import Seq

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
        node, log_score, match_type, base_in_lib, base_in_seq = x
        
        if node.i_seqset != prev_i_seqset and prev_i_seqset is not None:
            ideal_seq.append("*")
            real_seq.append("*")
        prev_i_seqset = node.i_seqset
            
        ideal_seq.append(base_in_lib)
        real_seq.append(base_in_seq)
    ideal_seq = "".join(ideal_seq)
    real_seq = "".join(real_seq)
    return ideal_seq, real_seq


def parse_node(node_str):
    assert node_str.startswith("Node(") and node_str.endswith(")")
    x = node_str[5:-1]
    x = x.split(",")
    x = [x.strip() for x in x if x.strip()]
    assert len(x) in [1, 4]
    return Node(*x)


def read_library(filename):
    "read_library(filename) -> list of (name, is_barcode, is_random, seqset)"
    import filelib

    # Format of file is:
    # "Name"  "Bar Code"  "Random Region"  "Alternate"  ["Alternate", ...]
    # <name>    <1/0>         <1/0>          <alt1>  <alt2>  <altN>
    # <name>    <1/0>         <1/0>          <alt1>  <alt2>  <altN>
    #
    # Columns are separated by tabs.  Each line contains a seqset.  A
    # seqset is represented as a list of sequences.
    #
    # Should have only 1 Bar Code.

    handle = filelib.read_cols(filename)
    header = handle.next()

    # Check the header format.
    assert len(header) >= 4, "invalid library format"
    assert header[0].upper() == "Name".upper()
    assert header[1].upper() == "Bar Code".upper()
    assert header[2].upper() == "Random Region".upper()
    for i in range(3, len(header)):
        assert header[i].upper() == "Alternate".upper()

    barcode_seen = False
    library = []  # list of (is_random, name, seqset)
    for i, cols in enumerate(handle):
        assert len(cols) >= 4
        name = cols[0]
        is_barcode = cols[1]
        is_random = cols[2]
        alternates = [x.strip() for x in cols[3:] if x.strip()]
        assert is_barcode in ["0", "1"]
        assert is_random in ["0", "1"]
        assert len(alternates) >= 1
        is_barcode, is_random = int(is_barcode), int(is_random)
        if is_barcode:
            assert not barcode_seen, "Can only have 1 barcode."
            barcode_seen = True
        x = SequenceSet(name, is_barcode, is_random, alternates)
        library.append(x)
    assert len(library)

    return library


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
        x = [handle.readline() for x in range(4)]
        lines = [x.strip() for x in x]
        if not lines[0]:
            break
        title, sequence, x, quality = lines
        assert x == "+"
        assert len(sequence) == len(quality)
        assert quality

        yield title, sequence, quality

