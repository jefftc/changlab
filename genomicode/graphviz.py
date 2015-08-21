"""

Uses the 1st quadrant in the Cartesian coordinate system. (0, 0, 0) is
at the lower left corner of the screen, extending positive to the
right and up.

TODO: fix this API

Functions:
make_graph    Creates an pygraphviz.AGraph object.
layout

"""

def make_graph(
    nodes, edges, node2attributes=None, edge2attributes=None,
    prog=None, subgraphs=None, directed=False, rank=False):
    # nodes is a list of the names of the nodes, given as strings.
    # edges is a list of tuples (<node_a>, <node_b>) that indicate
    # which nodes are connected.  node2attributes is a dictionary
    # where the key is the name of the node and the value is a
    # dictionary of attributes.  edge2attributes is a dictionary where
    # the key is a tuple in edges, and the value is a dictionary of
    # attributes.  prog is typically "neato" or "dot".  subgraph is
    # dictionary where the key is the name of the subgraph and value
    # is the list of nodes in that subgraph.  rank is to set
    # nodes in a subgraph to be same level,rank='same'
    #
    # G = make_graph(...)
    # G.draw(filename)
    # G.write(filename)
    
    # Node attributes:
    # style      filled
    # shape      box, circle, ellipse, point, triangle, diamond, octagon
    #            note, tab, folder
    # fillcolor  Color of background.  (style must be filled).
    # color      Color of outline, #FFFFFF
    #
    # Edge attributes:
    # style      dotted, bold
    # len        length
    # arrowhead
    import pygraphviz as pgv

    node2attributes = node2attributes or {}
    edge2attributes = edge2attributes or {}

    # Uses "neato" by default.
    prog = prog or "neato"
    subgraphs = subgraphs or {}

    forcelabels = True

    G = pgv.AGraph(dim=2, directed=directed)
    if forcelabels:
        G.graph_attr["forcelabels"] = "true"
    for node in nodes:
        attr = node2attributes.get(node, {})
        G.add_node(node, **attr)
    for i, j in edges:
        attr = edge2attributes.get((i, j), {})
        G.add_edge(i, j, **attr)
    for name in subgraphs:
        G.add_subgraph(subgraphs[name], name, rank=rank)
    G.layout(prog=prog)
    return G


def layout(nodes, edges, prog=None, subgraphs=None):
    # Return a list of the (x, y) coordinates, parallel to nodes.
    G = make_graph(nodes, edges, prog=prog, subgraphs=subgraphs)
    coords = []
    for node in nodes:
        n = G.get_node(node)
        x = n.attr["pos"]
        x, y = x.split(",")
        x, y = float(x), float(y)
        coords.append((x, y))
    return coords
