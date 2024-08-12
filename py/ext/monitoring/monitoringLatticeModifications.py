"""
Module. Includes a base function that will modify the accelerator lattice by inserting the monitoring accelerator nodes.
Based on setSC_General_AccNodes.
"""
# import the auxiliary classes
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the monitoring accelerator node constructor
from ext.monitoring.monitoringAccNodes import Monitoring_Base_AccNode

def setMonitoring_General_AccNodes(lattice, monitoring_ele_names,
        monitoring_length_min, monitoring_calculator,
        monitoring_NodeConstructor=None, insert_middle=False, monitoring_ele_types='any'):
    """
    It will put a set of a monitoring nodes into the lattice as child nodes of the first level accelerator nodes.
    The monitoring nodes will be inserted at the beginning of the part closest to the middle of the first level
    AccNode element. The function will return the array of monitoring nodes as a convenience for the user. This
    is a general function, and monitoring nodes will need specific information on what is actually measured.
    """

    monNodes_arr = [] # List of monitoring nodes that will be added

    # Identify the scattering calculator type
    monitoring_type = type(monitoring_calculator).__name__

    # Utility function to construct the monitoring node with given details
    def add_mon_node(lattNode, part_index):
        if monitoring_NodeConstructor == None:
            # Use the monitoring node constructor to make a lattice node
            monNode = Monitoring_Base_AccNode(monitoring_calculator, monitoring_type, 
                lattNode.getName()+":"+str(part_index)+":"+monitoring_type)
        else:
            monNode = monitoring_NodeConstructor(monitoring_calculator,
                    lattNode.getName()+":"+str(part_index)+":"+monitoring_type)

        monNodes_arr.append(monNode) # Append to the list of monitoring nodes
        # Append the node as a child node of the selected element
        lattNode.addChildNode(monNode,lattNode.BODY,part_index,lattNode.BEFORE)

    accNodes = lattice.getNodes() # Get all the nodes
    if(len(accNodes) == 0): return # Nothing to do if there aren't any nodes

    # First create an array of nodes to attach monitoring nodes to.
    # nodes_arr[(accNode, part_index, position, path_length)] 
    nodes_arr = [] # Empty array of nodes
    length_total = 0.
    running_path = 0.
    rest_length = 0.
    for accNode in accNodes: # Loop over all nodes
        nParts = accNode.getnParts() # Get the number of parts of each element
        for ip in range(nParts): # Loop over all parts
            part_length = accNode.getLength(ip) # Part length
            if(running_path > monitoring_length_min): # Append to the list of nodes to use 
                nodes_arr.append((accNode,ip,length_total,running_path))
                running_path = 0.
            running_path += part_length
            length_total += part_length
    if(len(nodes_arr) > 0): # Any remaining length?
        rest_length = length_total - nodes_arr[len(nodes_arr) - 1][2]
    else:
        rest_length = length_total
    #the first scattering node in the beginning of the lattice
    nodes_arr.insert(0,(accNodes[0],0,0.,rest_length))


    # Now we put the requested monitoring nodes as a children of accNodes
    ele_index_beg = ele_index_end = 0 # Keep track of the current element
    ele_name = nodes_arr[0][0].getName() # Current element name
    ele_type = nodes_arr[0][0].getType() # Current element type
    for inode in range(len(nodes_arr)):
        # Get details of the current and next node in the list
        (accNode, part_index, position, path_length) = nodes_arr[inode]
        if inode < len(nodes_arr)-1:
            (accNodeNext, part_indexNext, positionNext, path_lengthNext) = nodes_arr[inode+1]
        else: # this is the last node!
            positionNext = accNode.getLength(part_index)

        # Does the next node have a different name or is this the last element?
        if accNodeNext.getName() != ele_name or inode == len(nodes_arr)-1:
            ele_index_end = inode # This is the end of the element
            # Now that we have completely identified an element let's process
            # Check whether all elements are selected, or this one is selected
            # Also verify the type of the element
            if (monitoring_ele_names == 'all' or\
                (len(monitoring_ele_names) and ele_name in monitoring_ele_names)) and\
                (monitoring_ele_types == 'any' or\
                (len(monitoring_ele_types) and ele_type in monitoring_ele_types)):

                if insert_middle: # Do we need to insert in the middle
			        # Find the ideal mid-point position in the element
			        ele_part_pos_middle = 0.5*(nodes_arr[ele_index_beg][2] + positionNext)
                else:
                    ele_part_pos_middle = 0 # Make sure the insert position is at the start of the element

                # By default the monitor will be inserted at the beginning
                mon_insert_index = ele_index_beg
                mon_insert_pos = nodes_arr[ele_index_beg][2]
                # Interate over the part indices
                for ele_part_inode in range(ele_index_beg, ele_index_end+1):
                    # Find the first part which comes after the "mid-point".
                    # If ele_part_pos_middle is 0, then this sets the insertion at the beginning of the element.
                    if nodes_arr[ele_part_inode][2] >= ele_part_pos_middle:
                        mon_insert_index = ele_part_inode # Mark the index
                        mon_insert_pos = nodes_arr[ele_part_inode][2] # Mark the position
                        break # We found it!
                # print('Found {} domain = {},{} m, mon-pos = {} m'.format(ele_name,
                #      nodes_arr[ele_index_beg][2], positionNext, mon_insert_pos))
                # Add the monitoring node to the lattice
                add_mon_node(nodes_arr[mon_insert_index][0], nodes_arr[mon_insert_index][1])

        # We are done processing the current element, let's get to the next one
        if inode < len(nodes_arr)-1:
            ele_index_beg = inode+1 # Mark the position of the new element
            ele_name = accNodeNext.getName() # Save the name
            ele_type = accNodeNext.getType() # Save the type

    return monNodes_arr
