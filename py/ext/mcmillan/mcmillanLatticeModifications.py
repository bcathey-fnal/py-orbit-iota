"""
Module. Includes the function that will modify the accelerator lattice by placing McMillan nodes.
Based on setSC_General_AccNodes.
"""
# import the auxiliary classes
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the mcmillan accelerator node constructor
from ext.mcmillan.mcmillanAccNodes import McMillan_AccNode

def setMcMillan_AccNodes(lattice, mcmillan_ele_names, mcmillan_calculator):
    """
    Searches for the specific element names and replaces its default tracker with the McMillan tracker.
    Returns an array of mcmillan nodes as a convenience for the user.
    """

    mcmNodes_arr = [] # List of monitoring nodes that will be added

    # Utility function to construct the mcmillan node and insert it
    def add_mcm_node(lattNode):
        # First get important information about the node
        index = lattice.getNodeIndex(lattNode)
        nparts = lattNode.getnParts()
        nodelength = lattNode.getLength()
        ks = lattNode.getParam('B')
        # Use the mcmillan node constructor to make a lattice node
        mcmNode = McMillan_AccNode(mcmillan_calculator, nodelength, ks, lattNode.getName())
        mcmNode.setnParts(nparts)
        mcmNodes_arr.append(mcmNode) # Append to the list of mcmillan nodes
        #print('{}, index={}, l= {}, nparts={}'.format(mcmNode.__str__(), index, 
        #      mcmNode.getLength(), mcmNode.getnParts()))
        accNodes[index] =  mcmNode

    accNodes = lattice.getNodes() # Get all the nodes
    if(len(accNodes) == 0): return # Nothing to do if there aren't any nodes

    for accNode in accNodes: # Loop over all nodes
        ele_name = accNode.getName() # Get the name
        # Does the name match?
        if len(mcmillan_ele_names) and ele_name in mcmillan_ele_names:
            add_mcm_node(accNode)

    lattice.initialize()


    return mcmNodes_arr

