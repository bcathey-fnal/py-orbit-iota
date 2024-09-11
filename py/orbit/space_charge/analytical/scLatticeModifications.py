"""
Module. Includes functions that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""
import numpy as np # Linear interpolation of twiss functions

# import acc. nodes
from scAccNodes import SCanalyticalAccNode

# import the general SC lattice modification function 
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes

def setSCanalyticalAccNodes(lattice, sc_path_length_min,
							space_charge_calculator, twiss):
	"""
	Put a set of SCanalyticalAccNode nodes into the lattice as child nodes of
	the first level accelerator nodes.

	Parameters:
	- lattice: TEAPOT lattice in which space-charge nodes will be inserted.
	- sc_path_length_min: Minimum path length between consecutive space-charge
						  nodes.
	- space_charge_calculator: Object which will calculate and apply
							   space-charge kicks.
	- twiss: Numpy array with shape (Npts, 11) containing linear Twiss
			 parameter values. The columns are:
		0. Position in lattice (m)
		1. Phase advance of eigenmode a (2pi)
		2. Betatron amplitude of eigenmode a (m)
		3. Alfa function of eigenmode a
		4. Phase advance of eigenmode b (2pi)
		5. Betatron amplitude of eigenmode b (m)
		6. Alfa function of eigenmode b
		7. Horizontal dispersion (m)
		8. Derivative of horizontal dispersion
		9. Vertical dispersion (m)
		10. Derivative of vertical dispersion

	Returns:
	- List of inserted space-charge nodes.
	"""
	# Add the analytical space-charge nodes into the lattice
	scNodes_arr = setSC_General_AccNodes(lattice, sc_path_length_min,
	space_charge_calculator, SCanalyticalAccNode)
	s = 0.0 # Position of space-charge nodes
	for scNode in scNodes_arr:
		# First pack the relevant Twiss function valus into a dictionary
		lattice_functions = {}
		lattice_functions["betax"] = np.interp(s, twiss[:,0], twiss[:,2])
		lattice_functions["betay"] = np.interp(s, twiss[:,0], twiss[:,5])
		lattice_functions["etax"] = np.interp(s, twiss[:,0], twiss[:,7])
		lattice_functions["etay"] = np.interp(s, twiss[:,0], twiss[:,9])
		# TODO: INCLUDE A CLOSED ORBIT CALCULATION !!!
		lattice_functions["orbitx"] = 0.0
		lattice_functions["orbity"] = 0.0
		# Save them into the space-charge node
		scNode.set_lattice_functions(lattice_functions)
		scNode.setName(scNode.getName()+"FrozenSC") # Edit the node name
		s += scNode.getLength() # update the position
	lattice.initialize() # initialize the lattice
	return scNodes_arr