"""
Module. Includes a base function that will modify the accelerator lattice by inserting the scattering accelerator nodes.
Based on setSC_General_AccNodes.
"""
# import the auxiliary classes
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

def setScattering_General_AccNodes(lattice, scattering_ele_names, scattering_length_min, scattering_calculator, scattering_NodeConstructor):
    """
    It will put a set of a scattering nodes into the lattice as child nodes of the first level accelerator nodes.
    The scattering nodes will be inserted at the beginning of a particular part of the first level AccNode element.
    The distance between scattering nodes should be more than scattering_length_min. The function will return 
    the array of scattering nodes as a convenience for the user. This is a general function, and scattering nodes will need 
    specific information on what process is actually simulated.
    """
    accNodes = lattice.getNodes() # Get all the nodes
    if(len(accNodes) == 0): return # Nothing to do if there aren't any nodes

    # First create an array of nodes to attach scattering nodes to.
    # nodes_arr[(accNode, part_index, position, path_length)] 
    nodes_arr = [] # Empty array of nodes
    length_total = 0.
    running_path = 0.
    rest_length = 0.
    for accNode in accNodes: # Loop over all nodes
        nParts = accNode.getnParts() # Get the number of parts of each element
        for ip in range(nParts): # Loop over all parts
            part_length = accNode.getLength(ip) # Part length
            if(part_length > 1.0): # Check part length
                print "Warning! Node ",accNode.getName(), " has length ", part_length, "m. Scattering algorithms may be inaccurate."
            if(running_path > scattering_length_min): # Append to the list of nodes to use 
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


    # Now we put all scattering nodes as a children of accNodes
    scNodes_arr = []
    coolerLengths_arr = []
    ncoolers = 0
    current_cooler_elename = ''
    for inode in range(len(nodes_arr)-1):
        # Get details of the current node in the list
        (accNode, part_index, position, path_length) = nodes_arr[inode]
        # If the node is not in a list of names provided then skip
        if len(scattering_ele_names) > 0 and accNode.getName() not in scattering_ele_names:
            continue

        # Get the next element
        (accNodeNext, part_indexNext, positionNext, path_lengthNext) = nodes_arr[inode+1]
        
        # Keep track of total length of each cooler element
        if current_cooler_elename == accNode.getName():
            coolerLengths_arr[ncoolers-1] += path_lengthNext
        else:
            current_cooler_elename = accNode.getName()
            coolerLengths_arr.append(path_lengthNext)
            ncoolers += 1

        #print('LatticeModify: {}:{} l = {} m'.format(accNode.getName(),part_index,path_lengthNext))

        # Construct the scattering node with given details
        scNode = scattering_NodeConstructor(scattering_calculator, accNode.getName()+":"+str(part_index)+":")
        scNode.setLengthOfScattering(path_lengthNext) # set the length of the element
        scNodes_arr.append(scNode) # Append to the list of scattering nodes
        # Append the node as a child node of the selected element
        accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)
    # Set the last scattering node if no nodes are selected
    if len(scattering_ele_names) == 0:
        (accNode, part_index, position, path_length) = nodes_arr[len(nodes_arr)-1]
        scNode = scattering_NodeConstructor(scattering_calculator, accNode.getName()+":"+str(part_index)+":")
        scNode.setLengthOfScattering(rest_length)
        scNodes_arr.append(scNode)
        accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)

    return scNodes_arr, coolerLengths_arr

