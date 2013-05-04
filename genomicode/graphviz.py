"""

Uses the 1st quadrant in the Cartesian coordinate system. (0, 0, 0) is
at the lower left corner of the screen, extending positive to the
right and up.

Functions:
layout

"""

def layout(nodes, edges, prog=None, subgraphs=None):
    # nodes is a list of the names of the nodes, given as strings.
    # edges is a list of tuples (<node_a>, <node_b>) that indicate
    # which nodes are connected.  subgraph is dictionary where the key
    # is the name of the subgraph and value is the list of nodes in
    # that subgraph.  Return a list of the (x, y) coordinates,
    # parallel to nodes.
    import pygraphviz as pgv

    # Uses "neato" by default.
    prog = prog or "neato"
    subgraphs = subgraphs or {}

    G = pgv.AGraph(dim=2)
    for node in nodes:
        G.add_node(node)
    for i, j in edges:
        G.add_edge(i, j)
    for name in subgraphs:
        G.add_subgraph(subgraphs[name], name)
    G.layout(prog=prog)
    #print G.subgraphs()
    #print G.string()

    coords = []
    for node in nodes:
        n = G.get_node(node)
        x = n.attr["pos"]
        x, y = x.split(",")
        x, y = float(x), float(y)
        coords.append((x, y))
    return coords
