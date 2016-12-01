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

    # Header:
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
    if not data:
        return
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
