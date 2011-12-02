"""

read_gtr_file
read_atr_file
read_gtc_file
read_atc_file
read_kgg_file
read_kag_file

write_gtr_file
write_atr_file
write_gtc_file
write_atc_file
write_kgg_file
write_kag_file

parse_node
format_node

cut_dendrogram

"""

def read_gtr_file(filename):
    # Reads gtr or atr files.
    # Return a list of (left_id, right_id, dist).
    # The dist is the joining correlation between the left and right
    # nodes.
    # The index of each element in the list determines the node id.
    # node id = -(index+1)
    tree = []
    for i, line in enumerate(open(filename)):
        x = line.rstrip().split()
        node_str, left_str, right_str, dist = x
        x = map(parse_node, (node_str, left_str, right_str))
        node_id, left_id, right_id = x
        dist = float(dist)
        assert node_id == -(i+1)
        x = left_id, right_id, dist
        tree.append(x)
    return tree

def read_atr_file(filename):
    return read_gtr_file(filename)

def read_gtc_file(filename):
    # Format is:
    # GENE0X  <cluster>
    data = {}
    for line in open(filename):
        x = line.strip().split()
        assert len(x) == 2
        s, cluster = x[:2]
        cluster = int(cluster)
        if cluster == -1:
            cluster = None
        id = parse_node(s)
        data[id] = cluster
    return data

def read_atc_file(filename):
    return read_gtc_file(filename)

def read_kgg_file(filename):
    import jmath 

    # <ID> GROUP
    #
    # GROUP should be integer.
    handle = open(filename)
    x = handle.readline().strip().split()
    assert len(x) == 2
    assert x[1] == "GROUP"
    cluster = []
    for line in handle:
        x = line.rstrip("\r\n").split("\t")
        id, group = x
        group = jmath.safe_int(group)
        x = id, group
        cluster.append(x)
    return cluster

def read_kag_file(filename):
    return read_kgg_file(filename)

def _write_gtr_file_h(data, handle, item):
    import math
    # Calculate the width of the column.
    # NODE 10 X
    colwidth = 4 + int(math.ceil(math.log(len(data), 10))) + 1
    
    for i, x in enumerate(data):
        left_id, right_id, dist = x
        node_id = -i-1
        node_str = format_node(node_id, "NODE")
        left_str = format_node(left_id, item)
        right_str = format_node(right_id, item)
        print >>handle, "%-*s %-*s %-*s %.6f" % (
            colwidth, node_str, colwidth, left_str, colwidth, right_str, dist)

def write_gtr_file(data, handle):
    # Format:
    # <node_str> <left_str> <right_str> <distance>
    _write_gtr_file_h(data, handle, "GENE")

def write_atr_file(data, handle):
    _write_gtr_file_h(data, handle, "ARRY")

def _write_gtc_file_h(data, handle, item):
    # data is id -> cluster
    ids = sorted(data)
    for id in ids:
        s = format_node(id, item)
        cluster = data[id]
        if cluster is None:
            cluster = -1
        print >>handle, "%s %d" % (s, cluster)

def write_gtc_file(data, handle):
    # Format:
    # <node_str> <cluster>
    _write_gtc_file_h(data, handle, "GENE")

def write_atc_file(data, handle):
    _write_gtc_file_h(data, handle, "ARRY")

def write_kgg_file(data, handle):
    # Bug: Does not preserve original ID name.
    # Format (has header):
    # ID     GROUP
    # <id>   <cluster>
    #
    # The <id> might be a node string or the IDs provided in the PCL
    # file.  cluster seems to use whichever IDs are in the first
    # column of the PCL file.
    x = "ID", "GROUP"
    print >>handle, "\t".join(x)
    for id, cluster in data:
        x = id, cluster
        if cluster is None:
            cluster = "NA"
        print >>handle, "\t".join(map(str, x))

def write_kag_file(data, handle):
    write_kgg_file(data, handle)

def parse_node(s):
    assert (s.startswith("NODE") or
            s.startswith("GENE") or s.startswith("ARRY")), s
    assert s.endswith("X")
    id = int(s[4:-1])
    if s.startswith("NODE"):
        id = -id
    return id

def format_node(id, item):
    # item should be NODE, GENE, or ARRY.
    if id < 0:
        s = "NODE%dX" % -id
    else:
        s = "%s%dX" % (item, id)
    return s

def cut_dendrogram(tree, k):
    # Use a BFS to cut the tree into k clusters.  Return a dictionary
    # of id -> cluster.
    assert len(tree)
    
    tree_cluster = {}  # id -> cluster
    depth = 1          # The root is at depth 1.
    next_cluster = 0

    # Do a BFS through the tree, assigning clusters.
    stack = []         # node id, dist, cluster_num
    
    # Assume the root node is the last element of the tree and add it
    # to the stack.
    stack.append((-len(tree), tree[-1][2], None))
    while stack:
        # Sort the nodes from highest (closest to root) to lowest.
        # Small numbers are high nodes, big numbers are low nodes.
        stack = sorted(stack, cmp=lambda x, y: cmp(x[1], y[1]))

        # To do a BFS, get the highest node.
        id, dist, num = stack.pop(0)

        # If deep enough, and there's no cluster, create a new one.

        # If this node is not already in a cluster, and I'm deep
        # enough down the tree (for the number of clusters wanted), or
        # this node is a leaf (which is automatically deep enough),
        # then assign it to a cluster.
        if num is None and (depth >= k or id >= 0):
            num = next_cluster
            next_cluster += 1
            assert next_cluster <= k, "more clusters than items"

        # Assign the cluster.
        tree_cluster[id] = num

        if id >= 0:
            # If this is a leaf, then do nothing more.
            continue

        # Add the left and right nodes to the stack.
        left_id, right_id, dist = tree[-id-1]
        left_num = right_num = num
        left_dist = right_dist = 0
        if left_id < 0:
            # If this is an internal node, then set the distance.  The
            # default distance for leaves are 0, so they'll be
            # processed before internal nodes.
            left_dist = tree[-left_id-1][2]
        if right_id < 0:
            right_dist = tree[-right_id-1][2]
        stack.append((left_id, left_dist, left_num))
        stack.append((right_id, right_dist, right_num))

        # I've gone down one node, so increase the depth.
        depth += 1

    return tree_cluster

##     # The clustering reordered the items in the expression data set.
##     # Order cluster so that it matches the order in the expression
##     # data.
##     # item_order  item_ids      oindex   cindex
##     # GENE12X     200648_s_at     12       0
##     # GENE5X      200604_s_at      5       1
##     # GENE20X     200722_s_at     20       2
##     # [...]
##     #
##     # item2cluster is based on oindex.  Need to return the clusters in
##     # cindex order.

##     # Calculate the original indexes.
##     oindex = [_parse_gene_or_node(x) for x in item_order]

##     assert len(item2cluster) == len(tree)+1
##     cluster = [None] * (len(item2cluster))
##     for oi, n in item2cluster.iteritems():
##         ci = oindex.index(oi)
##         cluster[ci] = item_ids[oi], n

##     return cluster

