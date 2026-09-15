"""
Module. Puts the tune buffer nodes into a lattice at chosen positions.
"""
# import the tune buffer node
from ext.monitoring.tunebuffer.tbAccNodes import TuneBuffer_AccNode


def partStarts(lattice):
    """
    The position (m) at which every part of every node of the lattice starts,
    in order, and the node and part index of each.

    Returns: a list of (s, accNode, part_index).
    """
    starts = []
    s = 0.0
    for accNode in lattice.getNodes():
        for ip in range(accNode.getnParts()):
            starts.append((s, accNode, ip))
            s += accNode.getLength(ip)
    return starts


def setTuneBufferAccNodes(lattice, tunebuffer, positions, tolerance=1e-9):
    """
    Puts a TuneBuffer_AccNode at the start of the part of the lattice which
    begins at each of the positions (m), as child nodes before the part's
    body. The node at the smallest position is marked as the first of the
    ring. Every position must be the start of a part to within the
    tolerance; otherwise a ValueError is raised.

    Returns the list of nodes added, in order of position.
    """
    targets = sorted(float(s) for s in positions)
    nodes = []
    itarget = 0
    for (s, accNode, ip) in partStarts(lattice):
        while itarget < len(targets) and abs(targets[itarget] - s) <= tolerance:
            node = TuneBuffer_AccNode(tunebuffer,
                    "{}:{}:TuneBuffer".format(accNode.getName(), ip),
                    isfirst=(len(nodes) == 0))
            accNode.addChildNode(node, accNode.BODY, ip, accNode.BEFORE)
            nodes.append(node)
            itarget += 1
    if itarget < len(targets):
        raise ValueError("no part of the lattice starts at s = {} m, where a"
                         " tune buffer node was asked for".format(
                             targets[itarget]))
    lattice.initialize()
    return nodes
