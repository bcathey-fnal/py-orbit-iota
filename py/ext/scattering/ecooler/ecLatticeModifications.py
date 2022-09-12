"""
Module. Includes functions that will modify the accelerator lattice by inserting the ecooler accelerator nodes.
"""

# import SC acc. nodes
from ext.scattering.ecooler import ECooler_AccNode

# import general accelerator elements and lattice
#from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker


#import the general scattering lattice modification function 
from ext.scattering.scatteringLatticeModifications import setScattering_General_AccNodes

# Math functions
import math

def setECoolerAccNodes(lattice, cooler_ele_name, scattering_length_min, scattering_calculator):
    """
    Puts a set of ECooler_AccNodes into the lattice as child nodes of selected lattice elements
    """

    # Use the general scattering node creator function
    scNodes_arr, coolerLengths_arr = setScattering_General_AccNodes(lattice,
            cooler_ele_name, scattering_length_min, scattering_calculator, ECooler_AccNode)    
    dE, dtheta, dphi = scattering_calculator.getErrors() # Get the errors for the ecooler 

    # Initialize the name and count of number of coolers in the lattice
    ele_name = ''; ele_cnt = 0
    for scNode in scNodes_arr: # loop over all ecooler nodes
        # get the name and length of the part
        node_name = scNode.getName()
        node_elename, partnum, dummy = node_name.split(':')
        node_length = scNode.getLengthOfScattering()
        
        # Does this node belong to the last identified ecooler element
        if node_elename == ele_name:
            # Calculate the transverse offset based on the misalignment angles
            # dtheta and dphi
            ele_offset_x = ele_running_length*math.sin(dtheta)*math.cos(dphi)
            ele_offset_y = ele_running_length*math.sin(dtheta)*math.sin(dphi)
        else: # We're entering a new cooler
            ele_cnt += 1 # Increment the number of ecoolers
            ele_name = node_elename # Note the element name
            # Initialize the offset and running length
            ele_offset_x = ele_offset_y = ele_running_length = 0.0

        ele_seglen = node_length/math.cos(dtheta) # The default length
        # On the last part of the element, the part length will be smaller
        if ele_running_length+ele_seglen > coolerLengths_arr[ele_cnt-1]:
            ele_seglen = coolerLengths_arr[ele_cnt-1] - ele_running_length
        # But if the electron beam is so misaligned that the cooler segment ends
        # even before we iterate to the last part, we say it's too much!
        if ele_seglen <= 0.0:
            raise ValueError('Electron beam misalignment dtheta is too large!')
            
        ele_running_length += ele_seglen # Update the running length
        scNode.setName(node_name+"ECooler") # Update name
        scNode.setLengthOfScattering(ele_seglen) # Update segment length
        scNode.setoffset(ele_offset_x, ele_offset_y) # Set the offest
        #print(scNode) # Debug printing

    # initialize the lattice
    lattice.initialize()
    return scNodes_arr
		
