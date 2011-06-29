"""

Functions:
layout

"""

def layout(nodes, edges):
    import pygraphviz as pgv

    # Uses "neato" by default.
    prog = prog or "neato"

    G = pgv.AGraph(dim=2)
    for node in nodes:
        G.add_node(node)
    for i, j in edges:
        G.add_edge(i, j)
    G.layout(prog=prog)

    coords = []
    for node in nodes:
        n = G.get_node(node)
        x = n.attr["pos"]
        x, y = x.split(",")
        x, y = float(x), float(y)
        coords.append((x, y))
    return coords
