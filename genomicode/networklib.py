"""

Functions:
find_shortest_path
find_shortest_path_many

write_sif
write_noa

"""

class Network:
    def __init__(self, nodes, edges):
        # A node should be unique and can be any hashable Python
        # object.  edges is a list of (node1, node2).
        # Make sure no duplicates in nodes.
        seen = {}
        for n in nodes:
            assert n not in seen, "Duplicate node: %s" % n
            seen[n] = 1
        # Make sure all edges contain known nodes.
        for n1, n2 in edges:
            assert n1 in nodes
            assert n2 in nodes
        self.nodes = nodes[:]
        self.edges = edges[:]
        self._cache = None  # node -> list of neighbors
    def _cache_neighbors(self):
        cache = {}  # node -> node -> 1
        for n1, n2 in self.edges:
            if n1 not in cache:
                cache[n1] = {}
            if n2 not in cache:
                cache[n2] = {}
            cache[n1][n2] = 1
            cache[n2][n1] = 1
        for n in self.nodes:
            if n not in cache:
                cache[n] = {}
        cache2 = {}  # node -> list of neighbors
        for n in cache:
            cache2[n] = sorted(cache[n])
        self._cache = cache2
    def get_neighbors(self, node):
        assert node in self.nodes, "Unknown node: %s" % node
        if not self._cache:
            self._cache_neighbors()
        return self._cache[node]
    def num_neighbors(self, node):
        x = self.get_neighbors(node)
        return len(x)

def find_shortest_path(network, node1, node2):
    # Do a breadth-first search for the shortest path from node1 to
    # node2.  Return a list of all the paths that are shortest.  If
    # there is no path from node1 to node2, return an empty list.
    node2dist = {}
    node2dist[node1] = 0
    stack = [[node1]]   # list of proposed paths
    found_paths = []    # list of paths found
    iter = 0
    while stack:
        #print iter, len(stack), len(stack[0])
        path = stack.pop(0)
        # If path is complete, then add it to my list.
        if path[-1] == node2:  # found a good path
            #assert not found_paths or len(path) == len(found_paths[0])
            found_paths.append(path)
            continue
        # If this path is too long, then stop searching this path.
        if found_paths and len(path) >= len(found_paths[0]):
            continue
        # If this path is not complete, then search along the next node.
        for n in network.get_neighbors(path[-1]):
            if n in node2dist and node2dist[n] < len(path): # no cycles
                continue
            node2dist[n] = len(path)
            p = path + [n]
            stack.append(p)
        iter += 1
        #if iter > 10:
        #    import sys; sys.exit(0)
    # Bug checking: make sure each path is the same length.
    for path in found_paths:
        assert len(path) == len(found_paths[0])
    return found_paths

def _find_shortest_path_many_h(network, node_pairs):
    shortest_paths = {}
    for (n1, n2) in node_pairs:
        x = find_shortest_path(network, n1, n2)
        shortest_paths[(n1, n2)] = x
    return shortest_paths
            
def find_shortest_path_many(network, node_pairs, NUM_PROCS=None):
    # node_pairs is a list of (node1, node2).  Returns a dictionary of
    # (node1, node2) -> list of shortest paths.
    import multiprocessing

    NUM_PROCS = NUM_PROCS or 1
    assert NUM_PROCS >= 1
    BATCH_SIZE = 50

    # Split node pairs up so that each batch is at most BATCH_SIZE.
    batches = []
    while node_pairs:
        x = node_pairs[:BATCH_SIZE]
        node_pairs = node_pairs[BATCH_SIZE:]
        batches.append(x)

    if NUM_PROCS is not None:
        NUM_PROCS = min(NUM_PROCS, len(batches))

    # Distribute batches across NUM_PROCS 
    pool = multiprocessing.Pool(NUM_PROCS)
    results = []
    for batch in batches:
        if NUM_PROCS == 1:
            x = _find_shortest_path_many_h(network, batch)
        else:
            x = pool.apply_async(
                _find_shortest_path_many_h, (network, batch), {})
        results.append(x)
    pool.close()
    pool.join()
    if NUM_PROCS > 1:
        for i in range(len(results)):
            results[i] = results[i].get()

    shortest_paths = {}
    for x in results:
        shortest_paths.update(x)
    return shortest_paths

def write_sif(handle, network, edge2relationship):
    # network is a Network object.  edge2relationship is a dictionary
    # of (node1, node2) -> relationship, where relationship is a
    # string for the SIF format.

    seen = {}  # which nodes are seen
    for n1, n2 in network.edges:
        if (n1, n2) not in edge2relationship:
            n1, n2 = n2, n1
        assert (n1, n2) in edge2relationship
        seen[n1] = 1
        seen[n2] = 1
        rel = edge2relationship[(n1, n2)]
        x = n1, rel, n2
        print >>handle, "\t".join(map(str, x))
    # Print out all the single nodes.
    for n in network.nodes:
        if n in seen:
            continue
        print >>handle, n

def write_noa(handle, attribute, node2attribute):
    # attribute is the name of the attribute.  It should be a string
    # with no spaces.  node_attributes should be a dictionary of
    # node_name -> attribute_value.
    assert type(attribute) is type("")
    assert " " not in attribute
    print >>handle, attribute
    for node in sorted(node2attribute):
        print >>handle, "%s = %s" % (node, node2attribute[node])

