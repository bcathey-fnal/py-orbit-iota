"""
Module. Includes functions that will modify the accelerator lattice by inserting the tune detector nodes.
"""

# import monitoring acc. nodes
from ext.monitoring.tunedetector import TuneDetector_AccNode
#import the general Monitoring lattice modification function 
from ext.monitoring.monitoringLatticeModifications import setMonitoring_General_AccNodes

# import general accelerator elements and lattice
#from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker
from orbit.teapot import TEAPOT_Lattice, TEAPOT_MATRIX_Lattice # Lattice objs
from orbit.matrix_lattice import BaseMATRIX # For extracting matrices

# Linear algebra to handle Floquet matrices
import numpy as np

##############################################################################
# Miscellaneous Functions
# Convert the internal matrix format to numpy
def to_numpy_matrix(orbit_matrix):
    a = np.zeros((7,7)) # Declare an empty numpy matrix
    # Populate the matrix
    for i in range(0,7):
        for j in range(0,7):
            a[i,j] = orbit_matrix.get(i,j)
    return a

# Calculate Floquet transformation matrix
def get_4D_floquet_matrix(M_one_turn):
    R = np.zeros((4,4)) # Rotation matrix in normalized phase space
    M = M_one_turn[0:4,0:4] # 4D transfer matrix
    wM, vM = np.linalg.eig(M) # Get transverse eigenvalues
    # Fill in the rotation matrix
    R[0:2,0:2] = [[np.real(wM[0]), np.imag(wM[0])],
                  [np.imag(wM[1]), np.real(wM[1])]]
    R[2:4,2:4] = [[np.real(wM[2]), np.imag(wM[2])],
                  [np.imag(wM[3]), np.real(wM[3])]]
    # Get the eigen vectors for the rotation matrix.
    # Eigenvalues are exactly same by construction.
    wR, vR = np.linalg.eig(R)
    # Finally generate T such that M = inv(T) @ R @ T
    T = np.dot(vR, np.linalg.inv(vM))
    #Tinv = np.dot(vM, np.linalg.inv(vR))
    return np.real(T) # Clip the imaginary parts, which are ~0

##############################################################################
# Lattice modification function
def setTuneDetectorAccNodes(lattice, tunedetector_instance, bunch):
    """
    Puts a set of TuneDetector_AccNodes into the lattice as child nodes of
    selected lattice elements.
    """

    # First construct a matrix lattice for linear optics. Note that Matrix
    # nodes are only constructed for TEAPOT lattice elements.
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    # Store the one turn matrix
    M_one_turn = to_numpy_matrix(matrix_lattice.getRingMatrix())
    TFloquet_beg = get_4D_floquet_matrix(M_one_turn) # Get the Floquet matrix

    # Use the general monitoring node creator function to make all the
    # tune detector nodes and save the array. The Floquet matrices of
    # each node will be set later.
    tdNodes_arr = setMonitoring_General_AccNodes(lattice, 'all', 0.01,
            tunedetector_instance, TuneDetector_AccNode)

    tdNode_arr_cnt = 0
    M_cumulative = np.eye(7) # This will hold the product of the matrices
    for matrixNode in matrix_lattice.getNodes(): # loop over all nodes
        # Check whether this node is indeed a matrix node. Matrix nodes are
        # only constructed for TEAPOT lattice elements.
        if(isinstance(matrixNode,BaseMATRIX) == True): # Is this a matrix?
            # Get the matrix corresponding to the node
            node_matrix = to_numpy_matrix(matrixNode.getMatrix())
            # Do we still have td nodes to match?
            if tdNode_arr_cnt < len(tdNodes_arr):
                # Extract the parts of the node names
                matrixNode_name_parts = matrixNode.getName().split('_')
                tdNode_name_parts = tdNodes_arr[tdNode_arr_cnt].getName().split(':')
                # Match the base name and the part index
                if matrixNode_name_parts[0] == tdNode_name_parts[0] and \
                    matrixNode_name_parts[1] == tdNode_name_parts[1]:
                    #print('{} matched with {}'.format(matrixNode.getName(),
                    #    tdNodes_arr[tdNode_arr_cnt].getName()))
                    
                    # Extract the inverse of cumulative matrix
                    M_cumulative_inv = np.linalg.inv(M_cumulative[0:4,0:4])
                    wM, vM = np.linalg.eig(M_cumulative[0:4,0:4]) # Get transverse eigenvalues
                    R = np.zeros((4,4)) # Rotation matrix in normalized phase space
                    # Fill in the rotation matrix - WRONG CALCULATION!!!
                    R[0:2,0:2] = [[np.real(wM[0]), np.imag(wM[0])],
                                  [np.imag(wM[1]), np.real(wM[1])]]
                    R[2:4,2:4] = [[np.real(wM[2]), np.imag(wM[2])],
                                  [np.imag(wM[3]), np.real(wM[3])]]
                    # T' = R @ T @ inv(M_cum) - Pretty much cheating!
                    #tdNodes_arr[tdNode_arr_cnt].setTFloquet(np.dot(np.dot(
                    #    R, TFloquet_beg), M_cumulative_inv))
                    tdNode_arr_cnt += 1 # Now we move to the next tune detector node

            # Finally multiply the current matrix to keep track of the
            # cumulative transfer matrix
            M_cumulative = np.dot(node_matrix,M_cumulative)


    # Debug printing
    #for tdNode in tdNodes_arr: # loop over all tune detector nodes
    #    print(tdNode)

    # initialize the lattice
    lattice.initialize()
    return tdNodes_arr
		
