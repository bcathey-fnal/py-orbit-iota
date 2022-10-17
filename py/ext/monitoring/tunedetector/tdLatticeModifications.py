"""
Module. Includes functions that will modify the accelerator lattice by
inserting the tune detector nodes.
"""

# import monitoring acc. nodes
from ext.monitoring.tunedetector import TuneDetector_AccNode
# import the general Monitoring lattice modification function 
from ext.monitoring.monitoringLatticeModifications import setMonitoring_General_AccNodes
# Transfer matrices to convert to normal mode coords for a linear coupled lattice
from ext.linearoptics.Linear_Coupled_Lattice import Linear_Coupled_Lattice

# import general accelerator elements and lattice
from orbit.matrix_lattice import BaseMATRIX # For extracting matrices

###############################################################################
# Lattice modification function
def setTuneDetectorAccNodes(lattice, tunedetector_instance, bunch, verbosity):
    """
    Puts a set of TuneDetector_AccNodes into the lattice as child nodes of
    selected lattice elements.
    """

    # First construct a linear coupled, lattice for linear optics. Note that
    # Matrix nodes are only constructed for TEAPOT lattice elements.
    lcl = Linear_Coupled_Lattice(lattice, bunch)

    # Use the general monitoring node creator function to make all the
    # tune detector nodes and save the array. The transform matrices
    # for each node will be added later.
    tdNodes_arr = setMonitoring_General_AccNodes(lattice, 'all', 0.01,
            tunedetector_instance, TuneDetector_AccNode)

    tdNode_arr_cnt = 0
    tdNode_arr_len = len(tdNodes_arr)
    for matrixNode in lcl.getNodes(): # loop over all nodes
        # Check whether this node is indeed a matrix node. Matrix nodes are
        # only constructed for TEAPOT lattice elements.
        if(isinstance(matrixNode,BaseMATRIX) == True): # Is this a matrix?
            # Do we still have td nodes to match?
            if tdNode_arr_cnt < tdNode_arr_len:
                # Extract the parts of the node names
                matrixNode_name_parts = matrixNode.getName().split('_')
                tdNode_name_parts = tdNodes_arr[tdNode_arr_cnt].\
                                    getName().split(':')
                # Match the base name and the part index
                if matrixNode_name_parts[0] == tdNode_name_parts[0] and \
                    matrixNode_name_parts[1] == tdNode_name_parts[1]:
                    # Obtain the mode transform matrix for the lattice element
                    tdNodes_arr[tdNode_arr_cnt].setTmode(
                                matrixNode.getParam('Vinv'))
                    if verbosity > 1:
                        print('{} matched with {}'.format(
                            matrixNode.getName(),
                            tdNodes_arr[tdNode_arr_cnt].getName()))
                    # Now we move to the next tune detector node
                    tdNode_arr_cnt += 1

    # initialize the lattice
    lattice.initialize()
    # Print the number of nodes added
    if verbosity:
        print('{} TuneDetector nodes added to lattice {}'.format(
                tdNode_arr_cnt, lattice.getName()))
    return tdNodes_arr
		
